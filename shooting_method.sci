clc
clear
function [coefficients, polynomial] = interpolate(x, y)
    function y1 = vandermonde1(x1)
    n = length(x1)
    for i = 1:n
        y1(:, i) = x1^(i-1)
    end
    endfunction

    V = vandermonde1(x)
    coeffi = inv(V)*y'
    
    coefficients = coeffi
    polynomial = poly(coeffi', 'x', 'coeff')
    
endfunction

function y_dd = get_double_derivative(x, y, dy_dx)
    y_dd = -y
endfunction


function [x_list, y_list] = solve_ivp_forward(init_x, final_x, init_y, init_dy_dx)
    e = 0.01;
    n = (final_x/e) + 1 ;
    y_list = [];  
    a = get_double_derivative(init_x, init_y, init_dy_dx)
    whole_y = init_y;
    y_list($+1) = whole_y;  
    whole_dy = init_dy_dx;   
    half_y = init_y + (init_dy_dx * (e / 2));
    y_list($+1) = half_y;   
    half_dy = init_dy_dx + (a * (e / 2));   
    x = e;   
    while x <= final_x
        whole_y = whole_y + (half_dy * e);
        y_list($+1) = whole_y;
        
        a = get_double_derivative(x, half_y, half_dy);
        
        whole_dy = whole_dy + (a*e);
        
        
        half_y = half_y + (whole_dy*e);
        y_list($+1) = half_y;
        
        a = get_double_derivative(x, whole_y, whole_dy);
        
        half_dy = half_dy + (a*e);
        
        x = x + e;              
    end        
    x_list = linspace(0, final_x, length(y_list))  
    
    y_list = y_list'
endfunction

init_x = 0
final_x = %pi/2
init_y = 1
final_y = 1

a = -5 : 0.1 : 5;
F = [];

for i = -5 : 0.1 : 5
    [x1, y1] = solve_ivp_forward(init_x, final_x, init_y, i);
    F($+1) = y1(1, $) - final_y;    
end

F = F'
[coeffi, p] = interpolate(a, F)

exact_a = real(roots(a))(1, 1);

[x, y] = solve_ivp_forward(init_x, final_x, init_y, exact_a);
[xp, yp] = solve_ivp_forward(init_x, final_x, init_y, exact_a+0.5);
[xm, ym] = solve_ivp_forward(init_x, final_x, init_y, exact_a-0.5);

plot(x, y, 'r')
plot(xp, yp, 'k-')
plot(xm, ym, 'k-')
plot([final_x], [final_y], 'b.')
xlabel('x')
ylabel('y(x)')
h1 = legend(['best fit with a = ' + string(exact_a); 'with a = ' + string(exact_a + 0.5); 'with a = ' + string(exact_a -0.5)], 'in_upper_left')







