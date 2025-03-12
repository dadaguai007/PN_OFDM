clc;clear;
% 测试循环卷积矩阵


% 创建一个向量
vector = [1; 2; 3; 4];

% 创建一个对角矩阵，对角元素为 [5, 6, 7, 8]
diag_matrix = diag([5, 6, 7, 8]);

% 对角矩阵与向量相乘
result = diag_matrix * vector;

% 获取结果向量的长度
N = length(result);

% 初始化循环卷积矩阵
circulant_matrix = zeros(N, N);

% 构建循环卷积矩阵
for i = 1:N
    % 循环移位操作
    circulant_matrix(:, i) = circshift(result', i - 1);
end

% 显示原始向量、对角矩阵、相乘结果和循环卷积矩阵
disp('原始向量:');
disp(vector);

disp('对角矩阵:');
disp(diag_matrix);

disp('相乘结果:');
disp(result);

disp('循环卷积矩阵:');
disp(circulant_matrix);