clear all;
close all;

fig = figure;

% N = 100
T = [1: 16];
sequential_average = 0.000023267;
% sequential = [0.004426800  ];
openmp = [ 0.000036160  0.000235487  0.000327907  0.000085293  0.000372973  0.000547773  0.000188833  0.001139080  0.000394840  0.000449367  0.000417520  0.000401820  0.000482160  0.000443580  0.000467340  0.000537400 ];
subplot(2,2,1);  
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('100 x 100 Matrix');

% N = 200
subplot(5,2,2);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 0.000096460;
openmp = [ 0.000096487  0.000481413  0.000390467  0.000144633  0.000547853  0.000512660  0.001118667  0.002677233  0.001237327  0.000400673  0.000417900  0.000470287  0.000468640  0.000653087  0.000519313  0.000535113 ];
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('200 x 200 Matrix');

% N = 400
subplot(5,2,3);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = .000385900;
% sequential = [0.004426800  ];
openmp  = [0.000414813  0.000516633  0.000193467  0.000221693  0.000231807  0.000353533  0.000882140  0.001172500  0.000723953  0.000402280  0.000421047  0.000407187  0.000677433  0.000867933  0.000545760  0.000649980  ];
            

plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('400 x 400 Matrix');

% N = 800
subplot(5,2,4);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 0.001528827 ;
openmp  = [  0.001730347  0.001293467  0.000607887  0.000611200  0.000611547  0.000456247  0.000902653  0.002128653  0.001417107  0.000895600  0.000727227  0.000687173  0.000651553  0.000729587  0.000626207  0.000708787  ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('800 x 800 Matrix');




% N = 1600
subplot(5,2,5);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 0.006138173 ;
openmp  = [  0.006312853  0.003790627  0.002641560  0.002258813  0.002317713  0.002012767  0.002655273  0.003131153  0.002416727  0.002158080  0.002305947  0.001711813  0.001809873  0.001925080  0.002025220  0.002053280  ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('1,600 x 1,600 Matrix');

% N = 3200
subplot(5,2,6);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average =  0.036988533 ;
openmp  = [  0.033109260  0.019546300  0.013728733  0.011930827  0.012010933  0.010507820  0.009238420  0.011942620  0.009474920  0.009771433  0.008813600  0.009756633  0.008512553  0.009195400  0.008668080  0.010327220 ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('3,200 x 3,200 Matrix');


% N = 6400
subplot(5,2,7);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 0.141565740 ;
openmp  = [ 0.133694907  0.071922200  0.050018587  0.044008913  0.038722347  0.035935600  0.036055847  0.037481813  0.037422193  0.035794473  0.034771287  0.034942193  0.033091127  0.031369660  0.032244460  0.031418567  ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('6,400 x 6,400 Matrix');


% N = 12800
subplot(5,2,8);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 0.565845740 ;
openmp  = [ 0.560160493  0.303866927  0.218563120  0.187562533  0.163600167  0.151978620  0.136297200  0.128344113  0.144378780  0.129846960  0.131005213  0.128208733  0.129624273  0.129490267  0.125166213  0.122421420  ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('12,800 x 12,800 Matrix');


% N = 25600
subplot(5,2,9);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 2.106856647 ;
openmp  = [ 2.098523660  1.162218013  0.858291947  0.727438300  0.628460760  0.576268133  0.534719593  0.504219787  0.510848680  0.498216960  0.500837547  0.484737587  0.483719167  0.481168873  0.470412827  0.465384813   ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('25,600 x 25,600 Matrix');


% N = 51200
subplot(5,2,10);           
hold on;
grid on;
xlim([1 16]);
ylim([0 6]);
sequential_average = 6.619208547 ;
openmp  = [ 6.778169020  3.630253393  2.440990100  1.914836827  1.634250373  1.466635760  1.447087347  1.424355913  1.488720647  1.419839993  1.410212807  1.413176973  1.395706087  1.396149540  1.376563927  1.393046453  ];
            
plot(T, openmp.^-1 * sequential_average, 'LineWidth', 2);
title('51,200 x 51,200 Matrix');




handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup');
xlabel(handle,'Number of threads');
title(handle,'Speedup of parallelising method to count degrees of nodes in the graph');
handle.Title.Position(2) = 1.05
grid on;
    