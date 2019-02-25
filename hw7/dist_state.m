function Distn = dist_state(a_0, trans_P, T)
L = 30;
a = a_0 * trans_P ^T;

Distn_temp = zeros(L,L);
for i = 1:L
    for j = 1:L
        Distn_temp(i,j) = a(1, L *(i-1)+j);
    end
end

Distn = Distn_temp/(sum(sum(Distn_temp)));

end