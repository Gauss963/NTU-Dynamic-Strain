#!/usr/bin/env python3
import argparse
import math
import heapq
import numpy as np
import gmsh

# python3 check_tet_slivers.py ../Models/50mm-BS-PMMA.msh --topk 50 --rthr 1e-3 --hthr 1e-3

def build_node_indexer(node_tags: np.ndarray):
    node_tags = node_tags.astype(np.int64, copy=False)
    max_tag = int(node_tags.max())
    n = len(node_tags)
    if max_tag <= 10_000_000 and max_tag <= 5 * n:
        idx = np.full(max_tag + 1, -1, dtype=np.int64)
        idx[node_tags] = np.arange(n, dtype=np.int64)
        return ("dense", idx)
    d = {int(t): i for i, t in enumerate(node_tags)}
    return ("dict", d)


def map_tags_to_indices(mapper, tags_2d: np.ndarray):
    kind, obj = mapper
    if kind == "dense":
        return obj[tags_2d]
    out = np.empty_like(tags_2d, dtype=np.int64)
    it = np.nditer(tags_2d, flags=["multi_index"])
    while not it.finished:
        out[it.multi_index] = obj.get(int(it[0]), -1)
        it.iternext()
    return out


def tri_area(p0, p1, p2):
    # 0.5 * || (p1-p0) x (p2-p0) ||
    v1 = p1 - p0
    v2 = p2 - p0
    cr = np.cross(v1, v2)
    return 0.5 * np.linalg.norm(cr, axis=-1)


def write_pos_points(path: str, points_xyz: np.ndarray, values: np.ndarray, view_name="slivers"):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f'View "{view_name}" {{\n')
        for (x, y, z), v in zip(points_xyz, values):
            f.write(f"  SP({x},{y},{z}){{{v}}};\n")
        f.write("};\n")


def main():
    ap = argparse.ArgumentParser(description="Detect sliver tetra (tet4) by inradius/min altitude.")
    ap.add_argument("msh", help="Path to .msh file")
    ap.add_argument("--topk", type=int, default=50, help="Keep top-K worst tets")
    ap.add_argument("--chunk", type=int, default=200_000, help="Chunk size")
    ap.add_argument("--rthr", type=float, default=None, help="Count tets with inradius < rthr (mesh units)")
    ap.add_argument("--hthr", type=float, default=None, help="Count tets with min altitude < hthr (mesh units)")
    ap.add_argument("--pos_r", type=str, default="worst_inradius.pos", help="Output .pos for worst inradius")
    ap.add_argument("--pos_h", type=str, default="worst_altitude.pos", help="Output .pos for worst altitude")
    args = ap.parse_args()

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(args.msh)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags, dtype=np.int64)
    coords = np.asarray(node_coords, dtype=np.float64).reshape(-1, 3)
    mapper = build_node_indexer(node_tags)

    # get tet4 only: element type 4 in Gmsh
    types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(3)

    tet_tags = None
    tet_conn = None
    total_3d = 0
    for etype, etags, enodes in zip(types, elem_tags_list, elem_node_tags_list):
        props = gmsh.model.mesh.getElementProperties(etype)
        name = props[0]
        dim = props[1]
        if dim != 3:
            continue
        total_3d += len(etags)
        if etype == 4:  # Tetrahedron 4
            tet_tags = np.asarray(etags, dtype=np.int64)
            tet_conn = np.asarray(enodes, dtype=np.int64).reshape(-1, 4)

    if tet_tags is None:
        gmsh.finalize()
        raise RuntimeError("No Tetrahedron 4 elements (etype=4) found.")

    gmsh.finalize()

    ne = tet_tags.size
    print(f"[INFO] Nodes: {len(node_tags)}")
    print(f"[INFO] Total 3D elements: {total_3d}")
    print(f"[INFO] Tet4 count: {ne}")

    # heaps for worst (smallest) r and h
    heap_r = []  # store (-r, tag, bx,by,bz, r, hmin, V)
    heap_h = []  # store (-hmin, ...)

    global_min_r = float("inf")
    global_min_h = float("inf")
    global_min_V = float("inf")
    meta_r = None
    meta_h = None
    meta_V = None

    cnt_r = 0
    cnt_h = 0

    for start in range(0, ne, args.chunk):
        end = min(ne, start + args.chunk)
        conn = tet_conn[start:end]
        tags = tet_tags[start:end]

        idx = map_tags_to_indices(mapper, conn)
        if np.any(idx < 0):
            bad = np.sum(idx < 0)
            raise RuntimeError(f"Found {bad} unmapped node tags in tet connectivity.")

        p = coords[idx]  # (chunk,4,3)
        p0, p1, p2, p3 = p[:,0,:], p[:,1,:], p[:,2,:], p[:,3,:]

        # Volume V = |(p1-p0) dot ((p2-p0) x (p3-p0))| / 6
        v1 = p1 - p0
        v2 = p2 - p0
        v3 = p3 - p0
        V6 = np.abs(np.einsum("ij,ij->i", v1, np.cross(v2, v3)))
        V = V6 / 6.0

        # Face areas
        A012 = tri_area(p0, p1, p2)
        A013 = tri_area(p0, p1, p3)
        A023 = tri_area(p0, p2, p3)
        A123 = tri_area(p1, p2, p3)
        A_sum = A012 + A013 + A023 + A123

        # Inradius r = 3V / A_sum (avoid division by zero)
        r = np.where(A_sum > 0, 3.0 * V / A_sum, 0.0)

        # Min altitude: h_i = 3V / A_face_i ; min altitude = min_i h_i = 3V / max(A_face_i)
        A_max = np.maximum.reduce([A012, A013, A023, A123])
        hmin = np.where(A_max > 0, 3.0 * V / A_max, 0.0)

        # barycenter
        bc = (p0 + p1 + p2 + p3) / 4.0

        # global minima
        i_r = int(np.argmin(r))
        if float(r[i_r]) < global_min_r:
            global_min_r = float(r[i_r])
            meta_r = (int(tags[i_r]), float(bc[i_r,0]), float(bc[i_r,1]), float(bc[i_r,2]),
                      float(r[i_r]), float(hmin[i_r]), float(V[i_r]))

        i_h = int(np.argmin(hmin))
        if float(hmin[i_h]) < global_min_h:
            global_min_h = float(hmin[i_h])
            meta_h = (int(tags[i_h]), float(bc[i_h,0]), float(bc[i_h,1]), float(bc[i_h,2]),
                      float(r[i_h]), float(hmin[i_h]), float(V[i_h]))

        i_v = int(np.argmin(V))
        if float(V[i_v]) < global_min_V:
            global_min_V = float(V[i_v])
            meta_V = (int(tags[i_v]), float(bc[i_v,0]), float(bc[i_v,1]), float(bc[i_v,2]),
                      float(r[i_v]), float(hmin[i_v]), float(V[i_v]))

        # thresholds
        if args.rthr is not None:
            cnt_r += int(np.count_nonzero(r < args.rthr))
        if args.hthr is not None:
            cnt_h += int(np.count_nonzero(hmin < args.hthr))

        # update heaps (keep smallest r/hmin)
        # candidates from this chunk
        k = min(args.topk, end - start)

        cand_r = np.argpartition(r, k - 1)[:k]
        for ii in cand_r:
            rr = float(r[ii])
            item = (-rr, int(tags[ii]), float(bc[ii,0]), float(bc[ii,1]), float(bc[ii,2]),
                    rr, float(hmin[ii]), float(V[ii]))
            if len(heap_r) < args.topk:
                heapq.heappush(heap_r, item)
            else:
                worst_rr = -heap_r[0][0]  # largest r among kept
                if rr < worst_rr:
                    heapq.heapreplace(heap_r, item)

        cand_h = np.argpartition(hmin, k - 1)[:k]
        for ii in cand_h:
            hh = float(hmin[ii])
            item = (-hh, int(tags[ii]), float(bc[ii,0]), float(bc[ii,1]), float(bc[ii,2]),
                    float(r[ii]), hh, float(V[ii]))
            if len(heap_h) < args.topk:
                heapq.heappush(heap_h, item)
            else:
                worst_hh = -heap_h[0][0]
                if hh < worst_hh:
                    heapq.heapreplace(heap_h, item)

    # sort heaps ascending by r/h
    worst_r = sorted([(-it[0], *it[1:]) for it in heap_r], key=lambda x: x[0])  # (r, tag, bx,by,bz, r, hmin, V) note r duplicated
    worst_h = sorted([(-it[0], *it[1:]) for it in heap_h], key=lambda x: x[0])

    print("\n========== REPORT (Tet4) ==========")
    if meta_r:
        tag, bx, by, bz, rr, hh, VV = meta_r
        print(f"[MIN] inradius r = {rr:.6e} at tet tag {tag}, bc=({bx:.6e},{by:.6e},{bz:.6e}), hmin={hh:.6e}, V={VV:.6e}")
    if meta_h:
        tag, bx, by, bz, rr, hh, VV = meta_h
        print(f"[MIN] altitude hmin = {hh:.6e} at tet tag {tag}, bc=({bx:.6e},{by:.6e},{bz:.6e}), r={rr:.6e}, V={VV:.6e}")
    if meta_V:
        tag, bx, by, bz, rr, hh, VV = meta_V
        print(f"[MIN] volume V = {VV:.6e} at tet tag {tag}, bc=({bx:.6e},{by:.6e},{bz:.6e}), r={rr:.6e}, hmin={hh:.6e}")

    if args.rthr is not None:
        print(f"[COUNT] #tets with r < {args.rthr} : {cnt_r}")
    if args.hthr is not None:
        print(f"[COUNT] #tets with hmin < {args.hthr} : {cnt_h}")

    print("\n[TOP-K] Worst by inradius r (smallest first):")
    for i, row in enumerate(worst_r[:args.topk], 1):
        rr = row[0]; tag = row[1]; bx,by,bz = row[2],row[3],row[4]
        hmin = row[6]; V = row[7]
        print(f"  #{i:02d}: r={rr:.6e}, tag={tag}, bc=({bx:.6e},{by:.6e},{bz:.6e}), hmin={hmin:.6e}, V={V:.6e}")

    print("\n[TOP-K] Worst by min altitude hmin (smallest first):")
    for i, row in enumerate(worst_h[:args.topk], 1):
        hh = row[0]; tag = row[1]; bx,by,bz = row[2],row[3],row[4]
        rr = row[5]; V = row[7]
        print(f"  #{i:02d}: hmin={hh:.6e}, tag={tag}, bc=({bx:.6e},{by:.6e},{bz:.6e}), r={rr:.6e}, V={V:.6e}")

    # write .pos files
    pts_r = np.array([[row[2], row[3], row[4]] for row in worst_r], dtype=float)
    vals_r = np.array([row[0] for row in worst_r], dtype=float)
    write_pos_points(args.pos_r, pts_r, vals_r, view_name="worst_inradius_r")

    pts_h = np.array([[row[2], row[3], row[4]] for row in worst_h], dtype=float)
    vals_h = np.array([row[0] for row in worst_h], dtype=float)
    write_pos_points(args.pos_h, pts_h, vals_h, view_name="worst_altitude_hmin")

    print(f"\n[OUTPUT] Wrote: {args.pos_r} and {args.pos_h}")
    print("         Open them in Gmsh (File -> Open) to see where slivers are.")


if __name__ == "__main__":
    main()