close all
n = 20;
I = 100;
for delta = 1:n
    Jsum = 0;
    Esum = 0;
    for i = 1:I
        run LSequalizer.m
        Jsum = Jsum+ Jmin;
        Esum = Esum+err;
    end

    Javg(delta) = Jsum /I ;
    Eavg(delta) = Esum/I;
    if (Esum ==0)
        delta
    end
end
figure('Name','Javg'), plot(Javg(1:n),'b-')
xlabel('delta'), ylabel('avg(Jmin)'), grid on;
figure('Name', 'Eavg'), plot(1:n, Eavg(1:n),'x-')
xlabel('delta'), ylabel('avg(err)'),grid on;