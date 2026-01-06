import argparse
import math
import heapq
import numpy as np
import gmsh


def edge_pairs_for_primary_nodes(element_name: str, n_primary: int):
    """
    Return a list of (i,j) edge pairs using only primary (corner) nodes.
    Covers common 3D elements.
    """
    name = element_name.lower()

    # Tetra (4 primary nodes)
    if "tetra" in name and n_primary >= 4:
        return [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    # Pyramid (5 primary nodes)
    if "pyramid" in name and n_primary >= 5:
        return [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (0, 4), (1, 4), (2, 4), (3, 4),
        ]

    # Prism / Wedge (6 primary nodes)
    if ("prism" in name or "wedge" in name) and n_primary >= 6:
        return [
            (0, 1), (1, 2), (2, 0),
            (3, 4), (4, 5), (5, 3),
            (0, 3), (1, 4), (2, 5),
        ]

    # Hexahedron (8 primary nodes)
    if ("hexa" in name or "hex" in name) and n_primary >= 8:
        return [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]

    # Fallback: if unknown, use all pairs among primary nodes (can be expensive)
    # Limit to at most 8 nodes to avoid combinatorial explosion
    if n_primary <= 8:
        pairs = []
        for i in range(n_primary):
            for j in range(i + 1, n_primary):
                pairs.append((i, j))
        return pairs

    return []


def build_node_indexer(node_tags: np.ndarray):
    """
    Fast tag->index mapping. Uses dense array if tags are not too sparse, else dict.
    """
    node_tags = node_tags.astype(np.int64, copy=False)
    max_tag = int(node_tags.max())
    n = len(node_tags)

    # If tags are extremely sparse, dense mapping may be too big.
    # Heuristic: allow dense mapping up to ~10 million or 5x node count.
    if max_tag <= 10_000_000 and max_tag <= 5 * n:
        idx = np.full(max_tag + 1, -1, dtype=np.int64)
        idx[node_tags] = np.arange(n, dtype=np.int64)
        return ("dense", idx)

    # Fallback to dict
    d = {int(t): i for i, t in enumerate(node_tags)}
    return ("dict", d)


def map_tags_to_indices(mapper, tags_2d: np.ndarray):
    kind, obj = mapper
    if kind == "dense":
        return obj[tags_2d]
    # dict fallback (slower but safe)
    out = np.empty_like(tags_2d, dtype=np.int64)
    # vectorized dict lookup isn't available; do row-wise
    it = np.nditer(tags_2d, flags=["multi_index"])
    while not it.finished:
        out[it.multi_index] = obj.get(int(it[0]), -1)
        it.iternext()
    return out


def write_pos_points(path: str, points_xyz: np.ndarray, values: np.ndarray, view_name="small_edges"):
    """
    Write a simple Gmsh .pos with scalar points (SP).
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write(f'View "{view_name}" {{\n')
        for (x, y, z), v in zip(points_xyz, values):
            f.write(f"  SP({x},{y},{z}){{{v}}};\n")
        f.write("};\n")


def main():
    ap = argparse.ArgumentParser(description="Check min edge length in a Gmsh .msh (3D elements).")
    ap.add_argument("msh", help="Path to .msh file")
    ap.add_argument("--topk", type=int, default=30, help="Keep top-K smallest elements")
    ap.add_argument("--threshold", type=float, default=None,
                    help="Count elements with min edge length < threshold (mesh units)")
    ap.add_argument("--chunk", type=int, default=200_000, help="Chunk size for processing elements")
    ap.add_argument("--pos", type=str, default=None, help="Output .pos file to locate smallest elements")
    ap.add_argument("--sample", type=int, default=200_000,
                    help="Random sample size for approximate quantiles (set 0 to disable)")
    args = ap.parse_args()

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(args.msh)

    # Nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags, dtype=np.int64)
    coords = np.asarray(node_coords, dtype=np.float64).reshape(-1, 3)
    mapper = build_node_indexer(node_tags)

    print(f"[INFO] Nodes: {len(node_tags)}")

    # 3D elements only
    types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(3)

    global_min2 = float("inf")
    global_min_meta = None  # (etype, elem_tag, min_edge, barycenter)

    # top-K smallest using a max-heap on min_edge^2 (store negative to simulate max-heap)
    heap = []  # items: (-min2, etype, elem_tag, min_edge, bx, by, bz)

    below_thr_count = 0
    total_3d = 0

    # For approximate quantiles: reservoir-like random sample of min edges
    rng = np.random.default_rng(12345)
    sample_vals = []

    for etype, elem_tags, elem_node_tags in zip(types, elem_tags_list, elem_node_tags_list):
        props = gmsh.model.mesh.getElementProperties(etype)
        # props: (name, dim, order, numNodes, localNodeCoord, numPrimaryNodes)
        name = props[0]
        dim = props[1]
        num_nodes = props[3]
        n_primary = props[5] if len(props) > 5 and props[5] > 0 else num_nodes

        if dim != 3:
            continue

        edges = edge_pairs_for_primary_nodes(name, n_primary)
        if not edges:
            print(f"[WARN] Skip unsupported element type {etype} ({name}), primary nodes={n_primary}")
            continue

        elem_tags = np.asarray(elem_tags, dtype=np.int64)
        elem_node_tags = np.asarray(elem_node_tags, dtype=np.int64)

        ne = elem_tags.size
        total_3d += ne
        print(f"[INFO] Element type {etype} ({name}), count={ne}, num_nodes={num_nodes}, primary={n_primary}")

        conn = elem_node_tags.reshape(ne, num_nodes)[:, :n_primary]  # only primary nodes

        # processing chunks
        for start in range(0, ne, args.chunk):
            end = min(ne, start + args.chunk)
            conn_chunk = conn[start:end]
            tags_chunk = elem_tags[start:end]

            idx = map_tags_to_indices(mapper, conn_chunk)
            if np.any(idx < 0):
                bad = np.sum(idx < 0)
                raise RuntimeError(f"Found {bad} unmapped node tags in elements; mesh may be corrupted.")

            pts = coords[idx]  # (chunk, n_primary, 3)

            # min edge^2 per element
            min2 = np.full((end - start,), np.inf, dtype=np.float64)
            for (a, b) in edges:
                d = pts[:, a, :] - pts[:, b, :]
                l2 = np.einsum("ij,ij->i", d, d)
                np.minimum(min2, l2, out=min2)

            # global min
            local_argmin = int(np.argmin(min2))
            local_min2 = float(min2[local_argmin])
            if local_min2 < global_min2:
                global_min2 = local_min2
                etag = int(tags_chunk[local_argmin])
                bc = pts[local_argmin].mean(axis=0)
                global_min_meta = (int(etype), etag, math.sqrt(local_min2), float(bc[0]), float(bc[1]), float(bc[2]))

            # threshold counting
            if args.threshold is not None:
                thr2 = args.threshold * args.threshold
                below_thr_count += int(np.count_nonzero(min2 < thr2))

            # update top-K heap
            # take a small candidate set from this chunk
            k = min(args.topk, end - start)
            cand_idx = np.argpartition(min2, k - 1)[:k]
            for ii in cand_idx:
                l2 = float(min2[ii])
                l = math.sqrt(l2)
                bc = pts[ii].mean(axis=0)
                item = (-l2, int(etype), int(tags_chunk[ii]), l, float(bc[0]), float(bc[1]), float(bc[2]))
                if len(heap) < args.topk:
                    heapq.heappush(heap, item)
                else:
                    # heap[0] has the most negative? Actually min-heap on -l2 => smallest -l2 is most negative => largest l2
                    # We want to keep smallest l2, so replace if current l2 is smaller than worst in heap
                    worst_l2 = -heap[0][0]
                    if l2 < worst_l2:
                        heapq.heapreplace(heap, item)

            # sampling for quantiles
            if args.sample and args.sample > 0:
                # randomly pick up to a few values from this chunk
                take = min(2000, end - start)
                pick = rng.choice(end - start, size=take, replace=False)
                sample_vals.extend(np.sqrt(min2[pick]).tolist())
                # cap sample size
                if len(sample_vals) > args.sample:
                    # downsample to args.sample
                    keep = rng.choice(len(sample_vals), size=args.sample, replace=False)
                    sample_vals = [sample_vals[i] for i in keep]

    gmsh.finalize()

    if total_3d == 0:
        print("[ERROR] No 3D elements found.")
        return

    min_edge = math.sqrt(global_min2)
    print("\n========== REPORT ==========")
    print(f"[RESULT] Total 3D elements processed: {total_3d}")
    if global_min_meta is not None:
        etype, etag, lmin, bx, by, bz = global_min_meta
        print(f"[RESULT] Global min edge length = {lmin:.6e} (mesh units)")
        print(f"         At element tag = {etag}, element type = {etype}")
        print(f"         Barycenter = ({bx:.6e}, {by:.6e}, {bz:.6e})")
    else:
        print(f"[RESULT] Global min edge length = {min_edge:.6e} (mesh units)")

    if args.threshold is not None:
        print(f"[RESULT] Elements with min edge < {args.threshold} : {below_thr_count}")

    if sample_vals:
        sv = np.sort(np.asarray(sample_vals, dtype=np.float64))
        q = lambda p: sv[int(p * (len(sv) - 1))]
        print("[APPROX] Min/1%/5%/50%/95% (sampled) of min-edge:")
        print(f"         {q(0.0):.6e}, {q(0.01):.6e}, {q(0.05):.6e}, {q(0.50):.6e}, {q(0.95):.6e}")

    # topK sorted ascending by edge length
    top = sorted([(-it[0], *it[1:]) for it in heap], key=lambda x: x[0])  # (l2, etype, tag, l, bx, by, bz)
    print("\n[TOP-K] Smallest elements by min edge length:")
    for rank, (l2, etype, etag, l, bx, by, bz) in enumerate(top, 1):
        print(f"  #{rank:02d}: lmin={l:.6e}, etype={etype}, tag={etag}, bc=({bx:.6e},{by:.6e},{bz:.6e})")

    if args.pos:
        pts = np.array([[bx, by, bz] for (_, _, _, _, bx, by, bz) in top], dtype=float)
        vals = np.array([l for (_, _, _, l, _, _, _) in top], dtype=float)
        write_pos_points(args.pos, pts, vals, view_name="min_edge_topk")
        print(f"\n[OUTPUT] Wrote Gmsh view file: {args.pos}")
        print("         Open in Gmsh to locate the problematic region (Post-processing view).")


if __name__ == "__main__":
    main()