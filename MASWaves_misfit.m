
function  e = MASWaves_misfit(c_t,c_curve0)

Q = length(c_t); % Number of data points in theoretical and experimental dispersion curves.

temp = 0;
for q = 1:Q
    temp = temp + sqrt((c_curve0(q) - c_t(q))^2)/c_curve0(q);
end
e = 1/Q * temp * 100;

disp(['Misfit: e = ', num2str(e), '%'])

end
