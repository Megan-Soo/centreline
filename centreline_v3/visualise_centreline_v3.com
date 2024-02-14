$subject = '005'

gfx read node centreline_v3 region airway
gfx read elem template region airway

gfx read node /home/msoo935/Code/MRI-to-model_platform_v1/subjects/$subject/exfile/$subject-left region left_lung
gfx read node /home/msoo935/Code/MRI-to-model_platform_v1/subjects/$subject/exfile/$subject-right region right_lung

gfx mod g_e '/airway' node_points glyph sphere size "3*3*3" mat cyan label cmiss_number
gfx mod g_e '/airway' element_points mat cyan
gfx mod g_e '/airway' lines line_width 2 mat cyan

gfx cre win

gfx mod g_e '/left_lung' node_points material magenta
gfx mod g_e '/right_lung' node_points mat magenta

gfx edi sce
