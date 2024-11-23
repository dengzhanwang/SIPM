% 创建一个新的图形窗口并设置大小
figure('Position', [100, 100, 1200, 400]);  % 设置图形窗口的位置和大小

% 循环从1到5，为每个数字创建一个图表
for i = 1:5
    % 创建子图，指定边距保留空间
    subplot('Position', [0.15+(i-1)*0.15 0.1 0.13 0.8]);  % 调整每个子图的位置和大小
    
    % 生成数据并绘制柱状图，使用多个数据系列以便创建图例
    data = rand(5, 1) * 10;  % 随机生成一些数据
    hold on;  % 保持当前图形，允许在上面继续绘图
    for j = 1:length(data)
        bar(j, data(j), 'DisplayName', ['Label ' num2str(j)]);  % 为每个条形赋予图例标签
    end
    hold off;  % 释放图形锁定
    
    % 标注数字
    ylim([0, max(data)*1.2]);  % 调整y轴范围，确保标签可见
    text(3, max(data)*1.1, num2str(i), 'FontSize', 14, 'HorizontalAlignment', 'center');
    
    % 设置标题（可选）
    title(['Graph ' num2str(i)]);
end

% 在图形窗口中添加一个共同的图例
legend('show', 'Location', 'best');  % 将图例放在子图外的右侧
saveas(gcf,'figure1.png')