for i = 1:30
    Opti_Record(i,:) = A_Opti_Record{i};
end

average = sum(abs(Opti_Record(:,18)+0.398*ones(30,1)))/30;
