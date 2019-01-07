tic

for jjj = 1:10
    result = solve_triangles_and_point0(gl, test, 0.6 , 1.3, 0.91, 1, 6000, 0.5, chord_su_b, sigma, tmax_su_c);
end
t = toc;
disp(['Time for iteration = ', num2str(t)]);


%disp([result.points.b])]
