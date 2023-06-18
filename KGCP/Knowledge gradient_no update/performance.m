for i = 1:30
    a = isempty(A_Opti_Record{i});
    if (a ==0)
        Opti_Record(i,:) = A_Opti_Record{i};
    end
end

average50 = sum(abs(Opti_Record(:,51)-1.0316*ones(30,1)))/30;
average100 = sum(abs(Opti_Record(:,101)-1.0316*ones(30,1)))/30;
average200 = sum(abs(Opti_Record(:,181)-1.0316*ones(30,1)))/30;