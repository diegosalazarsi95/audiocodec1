real_float=real/32768;
construido_float=construido/32768;
construido2_float = construido2/256;
promedio = 0;
for x = 1:64
    disp(x);
    error(x)= abs((real_float(x) - construido_float(x))/real_float(x))*100;
    disp(error(x));
    promedio = promedio +error(x);
end
   B=sort(error);
t=0:125e-6:7.875e-3;


figure()
plot(t,real_float,'-o',t,construido_float,'-x',t,construido2_float,'-x','LineWidth',1.2)
xlabel('Tiempo (s)')
ylabel('Magnitud')
title('Comparación de señal original y reconstruida')
legend('Audio original','Audio reconstruido', 'Audio reconstruid2')

%%Fs =8000;
%%audiowrite('real.wav',real_float,Fs);
%%audiowrite('construido.wav',construido_float,Fs);