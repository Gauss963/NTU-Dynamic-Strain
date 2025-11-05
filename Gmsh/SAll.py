import gmsh
import sys

def main():

    PMMA_THICKNESSES = [50, 100]
    mesh_size = 5
    # mesh_size = 20

    for PMMA_thickness in PMMA_THICKNESSES:

        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1) # 1 = Delaunay
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        # SubdivisionAlgorithm = 1 是 Transfinite 的必要條件
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) 
        gmsh.option.setNumber("Mesh.Binary", 0)

        gmsh.model.add("ContactModel")

        # 步驟 1：定義 3 個塊的幾何
        # 塊 1: moving-block
        blk1_origin = (0, 0, 0)
        blk1_dims = (200, 500, PMMA_thickness)
        
        # 塊 2: stationary-block (底部, 匹配部分)
        blk2_bottom_origin = (200, 0, 0)
        blk2_bottom_dims = (145, 500, PMMA_thickness) # Y=500, 匹配 blk1

        # 塊 3: stationary-block (頂部, 不匹配部分)
        blk2_top_origin = (200, 500, 0) # Y 從 500 開始
        blk2_top_dims = (145, 50, PMMA_thickness)  # Y=50 (550 - 500)

        # 步驟 2：僅創建幾何體
        blk1_tag = gmsh.model.occ.addBox(
            blk1_origin[0], blk1_origin[1], blk1_origin[2],
            blk1_dims[0], blk1_dims[1], blk1_dims[2]
        )
        
        blk2_bottom_tag = gmsh.model.occ.addBox(
            blk2_bottom_origin[0], blk2_bottom_origin[1], blk2_bottom_origin[2],
            blk2_bottom_dims[0], blk2_bottom_dims[1], blk2_bottom_dims[2]
        )

        blk2_top_tag = gmsh.model.occ.addBox(
            blk2_top_origin[0], blk2_top_origin[1], blk2_top_origin[2],
            blk2_top_dims[0], blk2_top_dims[1], blk2_top_dims[2]
        )

        # 步驟 3：Fragment (縫合) 所有塊
        # 這是實現共形網格的關鍵
        all_vols = [(3, blk1_tag), (3, blk2_bottom_tag), (3, blk2_top_tag)]
        # 將所有 3 個卷與彼此進行 fragment
        ov, ovm = gmsh.model.occ.fragment(all_vols, [])
        gmsh.model.occ.synchronize()

        # 步驟 4：獲取 Fragment 之後的 *新* 體積標籤
        # ovm[0] 是 blk1_tag 的映射 -> new_blk1_tag
        # ovm[1] 是 blk2_bottom_tag 的映射 -> new_blk2_bottom_tag
        # ovm[2] 是 blk2_top_tag 的映射 -> new_blk2_top_tag
        
        new_blk1_tag = ovm[0][0][1]
        new_blk2_bottom_tag = ovm[1][0][1]
        new_blk2_top_tag = ovm[2][0][1]
        
        new_vols = [new_blk1_tag, new_blk2_bottom_tag, new_blk2_top_tag]

        # 步驟 5：在 *新* 體積上設置 Transfinite (Hex 網格)
        all_faces_tags = set()
        for vol_tag in new_vols:
            # 設置體積為 Transfinite
            gmsh.model.mesh.setTransfiniteVolume(vol_tag)
            gmsh.model.mesh.setRecombine(3, vol_tag)
            
            # 獲取此體積的所有面
            faces = gmsh.model.getBoundary([(3, vol_tag)], oriented=False)
            for dim, tag in faces:
                all_faces_tags.add(tag)

        # 設置所有面 (包括內部) 為 Transfinite
        for tag in all_faces_tags:
            gmsh.model.mesh.setTransfiniteSurface(tag)
            gmsh.model.mesh.setRecombine(2, tag)

        # 設置全局網格大小
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        
        gmsh.model.occ.synchronize()

        # 步驟 6：手動定義物理群組
        # 這是最穩健的方法，使用 BoundingBox 來查找實體
        
        # 6a. 體積群組
        gmsh.model.addPhysicalGroup(3, [new_blk1_tag], 11)
        gmsh.model.setPhysicalName(3, 11, "moving-block")
        
        gmsh.model.addPhysicalGroup(3, [new_blk2_bottom_tag, new_blk2_top_tag], 21)
        gmsh.model.setPhysicalName(3, 21, "stationary-block")

        # 6b. CZM 介面 (關鍵！)
        # 這個面位於 x=200, y=[0, 500]
        # 我們需要為 *同一個* 幾何面定義 *兩* 個物理群組
        tol = 1e-3
        czm_face = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 0 - tol, 0 - tol,
            200 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if not czm_face:
            print("錯誤：未找到 CZM 介面")
        else:
            czm_face_tag = czm_face[0][1] # 獲取 (dim, tag) 中的 tag
            
            # 這是你的 CZM 介面
            # 你的原始代碼使用 15 和 24
            slave_tag = 15 # tag_prefix * 10 + 5 (moving-block-back)
            master_tag = 24 # tag_prefix * 10 + 4 (stationary-block-front)
            shear_tag = 55
            
            gmsh.model.addPhysicalGroup(2, [czm_face_tag], slave_tag)
            gmsh.model.setPhysicalName(2, slave_tag, "friction_slave")
            
            gmsh.model.addPhysicalGroup(2, [czm_face_tag], master_tag)
            gmsh.model.setPhysicalName(2, master_tag, "friction_master")

        # 6c. (可選) 其他外部邊界
        # 為了簡潔起見，我只添加了 CZM 面。
        # 你可以使用類似的 getEntitiesInBoundingBox 邏輯來標記
        # 'top', 'bottom', 'left', 'right' 等。
        
        # 例如: 'moving-block-front' (x=0)
        front_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            0 + tol, 500 + tol, PMMA_thickness + tol,
            dim=2
        )
        if front_face:
            gmsh.model.addPhysicalGroup(2, [front_face[0][1]], 14) # tag 14
            gmsh.model.setPhysicalName(2, 14, "moving-block-front")
        left_face = gmsh.model.getEntitiesInBoundingBox(
            0 - tol, 0 - tol, 0 - tol,
            200 + tol, 0 + tol, PMMA_thickness + tol,
            dim=2
        )
        if left_face:
            gmsh.model.addPhysicalGroup(2, [left_face[0][1]], 13) # tag 13
            gmsh.model.setPhysicalName(2, 13, "moving-block-left")
            
            
        back_faces = gmsh.model.getEntitiesInBoundingBox(
            345 - tol, 0 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if back_faces:
            back_face_tags = [tag for dim, tag in back_faces]
            gmsh.model.addPhysicalGroup(2, back_face_tags, 25)  # tag 25
            gmsh.model.setPhysicalName(2, 25, "stationary-block-back")
        else:
            print("⚠️ 無法找到 stationary-block-back")
        
        right_faces = gmsh.model.getEntitiesInBoundingBox(
            200 - tol, 550 - tol, 0 - tol,
            345 + tol, 550 + tol, PMMA_thickness + tol,
            dim=2
        )
        if right_faces:
            right_face_tags = [tag for dim, tag in right_faces]
            gmsh.model.addPhysicalGroup(2, right_face_tags, 26)  # tag 26
            gmsh.model.setPhysicalName(2, 26, "stationary-block-right")
        else:
            print("⚠️ 無法找到 stationary-block-right")
    
        
    
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA-CZM.msh")
        gmsh.write(f"../Models/{PMMA_thickness}mm-PMMA.brep")
        print(f"成功生成 {PMMA_thickness}mm-PMMA-CZM.msh")
        print("CZM 介面 (Slave/Master) 已在共形網格上創建。")
        
        gmsh.fltk.run()
        gmsh.finalize()

if __name__ == "__main__":
    main()