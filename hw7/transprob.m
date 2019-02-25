function prob = transprob(i, j, p) % transition probability to state w1 -1, w1, w1 + 1, where current state is (w1,w2)
     d = demand(p(i,j), p(j,i)); % prob of firm 1 makes a sale
     L = 30;
     if i == 1
         prob = [0, d * depreciation(i) + (1 - d), d * (1 - depreciation(i))];
     elseif i == L
         prob = [(1-d) * depreciation(i), (1-d) * (1 - depreciation(i)) + d * depreciation(i), 0];
     else
         prob = [(1-d) * depreciation(i), (1-d) * (1 - depreciation(i)) + d * depreciation(i), d * (1 - depreciation(i))];
     end
end