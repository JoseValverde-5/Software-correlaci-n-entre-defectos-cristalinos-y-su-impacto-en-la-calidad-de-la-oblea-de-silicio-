A = imread('Densidadmodificado.jpg');
B =  imread('Tiempodevida.jpg');
figure
imshowpair(A,B,'blend','Scaling','joint'); %Muestra las imagenes superpuestas
figure
imshowpair(A,B,'montage');%Comparacion en paralelo de las gráficas
beep;
