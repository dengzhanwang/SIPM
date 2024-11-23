function matchedIndex = matchClusterLabels(index1, index2)
    % 获取 index1 和 index2 的唯一标签
    labels1 = unique(index1,'stable');
    labels2 = unique(index2,'stable');
    
    % 初始化标签映射
    labelMapping = zeros(1, length(labels1));
    
    % 构建标签映射
    for i = 1:length(labels1)
        % 找到 index1 中的每个标签在 index2 中的对应标签
        labelMapping(i) = find(labels2 == labels1(i));
    end
    
    % 使用标签映射修改 index2 中的标签
    matchedIndex = index2;
    for i = 1:length(labels1)
        matchedIndex(index2 == labels2(i)) = labels1(i);
    end
end
