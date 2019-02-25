function UI = indicator(u)
U=abs(u);
UI=U;
UI(U>1)=0;
UI(U<=1)=1;

        
