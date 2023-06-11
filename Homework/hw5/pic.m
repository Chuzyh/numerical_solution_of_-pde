% 待绘制的 mu 值
mus = [0.8, 1.6, 2, 2.4];

for i = 1:length(mus)
    mu = mus(i);
    h = 1/50;
    p = 0:1:50;

    z = 1i * sin(2*pi * p*h) + cos(2*pi * p*h) - 1;
    plot(real(z), imag(z)); 
    hold on

    n = 50;
    A = zeros(n);
    for j = 1:n-2
        A(j,j) = 3;
        A(j + 1, j) = -4;
        A(j + 2, j) = 1;
    end
    A(n-1,n-1) = 3;
    A(n, n-1) = -4;
    A(n,n) = 3;
    A(1, n) = -4;
    A(1, n - 1) = 1;
    A(2, n) = 1;

    B = zeros(n);
    for j = 1:n-2
        B(j,j) = 1;
        B(j + 1, j) = -2;
        B(j + 2, j) = 1;
    end
    B(n-1,n-1) = 1;
    B(n, n-1) = -2;
    B(n,n) = 1;
    B(1, n) = -2;
    B(1, n - 1) = 1;
    B(2, n) = 1;

    v_a = eig(A);
    v_b = eig(B);
    z_p = -mu/2.*v_a + mu^2/2 .* v_b;
    sz = 10; 
    scatter(real(z_p), imag(z_p),sz,'filled'); 
    axis equal;
    title(['mu = ', num2str(mu)]);
    % 保存图像
    saveas(gcf, ['mu_', num2str(mu), '.png']);
    % 清除图像
    clf;
end