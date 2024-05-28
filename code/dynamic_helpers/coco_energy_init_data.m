function data = coco_energy_init_data(fdata, oid)
data.coll_seg = fdata.coll_seg;
data.xbp_idx = data.coll_seg.maps.xbp_idx;
data.T_idx = data.xbp_idx(end) + 1;
data.oid = oid;
data = coco_func_data(data);
end