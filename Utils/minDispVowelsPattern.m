function pattern = minDispVowelsPattern(vowels)

    diff = zeros(size(vowels, 1), size(vowels, 1));
    for i=1:size(vowels,1)-1
        for j=i+1:size(vowels,1)
            diff(i,j) = sqrt(sum((vowels(i,:)-vowels(j,:)).^2));
            diff(j,i) = diff(i,j);
        end
    end
    pattern = diff;
end