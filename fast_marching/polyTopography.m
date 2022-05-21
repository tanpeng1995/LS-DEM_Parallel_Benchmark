fid = fopen('~/Desktop/MPI_code/Impulse-based/levelset_topography_zoomin.dat', 'r');
lset = fscanf(fid, '%f');
disp("finish reading topography morphological file.");
lset = reshape(lset, [2000, 640, 200]);
lset = permute(lset, [2,1,3]);
[faces, verts] = isosurface(lset, 0);
disp("finish finding isosurface: zero level-set.");
[faces, verts] = reducepatch(faces, verts, 0.5);
TR = triangulation(faces, verts);
VN = vertexNormal(TR);
disp("finish triangulation.");
outfilename = '~/Desktop/MPI_code/Impulse-based/poly_topography_zoomin.dat';
outfile = fopen(outfilename, 'w');
fprintf(outfile, '%d\n', length(verts));
for v = 1:length(verts)
    fprintf(outfile,'%.3f %.3f %.3f\n', verts(v,1), verts(v,2), verts(v,3)); 
end
disp("finish writing vertices.");
fprintf(outfile, '%d\n', length(faces));
for f = 1:length(faces)
    fprintf(outfile,'%d %d %d\n', faces(f,1)-1, faces(f,2)-1, faces(f,3)-1); 
end
disp("finish writing faces.")
for v = 1:length(verts)
    fprintf(outfile, '%.3f %.3f %.3f\n', VN(v,1), VN(v,2), VN(v,3));
end
disp("finish writing vertex normals.")
fclose(fid);