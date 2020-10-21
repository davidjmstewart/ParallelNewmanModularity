clear all;
close all;

fig = figure;
matrix_sizes = [128, 256, 512, 1024, 2048, 4096, 8192, 16384];
% matrix_sizes = [1000000, 2000000, 3000000, 4000000];

parallel_timings = [
0.000018507  0.001021027  0.000678180  0.000921913  0.000737740  0.001464787  0.002827580  0.003398427;
0.000078073  0.001175627  0.000612927  0.000917800  0.001285820  0.002052087  0.003384833  0.003656487;
0.000372093  0.000895620  0.001040127  0.001447967  0.001542840  0.002251007  0.003703087  0.004495967;
0.001549093  0.001199733  0.001373393  0.001266227  0.001597087  0.002211240  0.002623260  0.003759127;
0.005406727  0.003373927  0.003193407  0.003181473  0.003062553  0.003217807  0.004143987  0.005191253;
0.022234473  0.012848833  0.009228453  0.007413647  0.006615660  0.006869607  0.008766620  0.009343153;
0.086336613  0.045515567  0.032609907  0.027009893  0.025305107  0.023423133  0.023199540  0.024284413;
0.334162927  0.174504847  0.127825160  0.106875500  0.095626047  0.092990533  0.090253667  0.092837313; ];

python_sequential_averages = [
1.5890000000001737e-05, 0.00032615999999999754, 0.0007261499999999977, 0.0029500099999999916, 0.011393629999999978, 0.04766399000000003, 0.18549439999999962, 0.7896764600000011
];

T = [1: 8];

sequential_averages = [
0.000018567 
0.000090707 
0.000598660 
0.001977800 
0.006093933 
0.023465120 
0.086270107 
0.340525380 
]


for i = 1:size(parallel_timings, 1)
    subplot(4,2,i); 
    hold on;
    grid on;
grid minor;
%     xlim([1 16]);
%     ylim([0 6]);
    hold on;
    c_par = plot(T, parallel_timings(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);
    p_par = plot(T, parallel_timings(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);
    hold off;
    title(sprintf('%d element Vector',matrix_sizes(i)),'interpreter','latex','fontsize',12);
    set(gca,'TickLabelInterpreter','latex')

end

hL = legend([c_par,p_par],{'Speedup (relative to sequential C)','Speedup (relative to sequential Python)'},'interpreter','latex','fontsize',12);
handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.XLabel.Interpreter = 'latex';


handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup','interpreter','latex','fontsize',16);
xlabel(handle,'Number of threads','interpreter','latex','fontsize',16);
title(handle,'\textbf{Speedup of parallelising Matrix-Vector multiplication with OpenMP threads}','interpreter','latex','fontsize',12);
handle.Title.Position(2) = 1.05





fig = figure;
% matrix_sizes = [1000000, 2000000, 3000000, 4000000];

parallel_timings_simd = [
0.000103773  0.000832113  0.000845967  0.000542320  0.000819707  0.001177120  0.002609780  0.003318793;
0.000036107  0.000985560  0.000550040  0.000759460  0.000961893  0.001192100  0.001726460  0.003313367;
0.000142293  0.000614207  0.000967287  0.000890600  0.000740060  0.001457260  0.002874053  0.002970927;
0.000559207  0.000560133  0.000936607  0.001128687  0.001466100  0.001778240  0.003020340  0.003458947;
0.002810820  0.002225607  0.001941380  0.001715787  0.001862553  0.003023380  0.004615793  0.004897660;
0.010836700  0.007781313  0.006945493  0.006026313  0.006113440  0.005878860  0.006390147  0.007121860;
0.036621093  0.025393027  0.021168053  0.019912653  0.019434240  0.019027033  0.019665740  0.020057827;
0.143174393  0.102118080  0.087503013  0.078300853  0.074452900  0.071914560  0.073481580  0.077631493;]

T = [1: 8];


parallel_timings_with_cache_block_64 = [
    
0.000012160  0.001086327  0.001303287  0.000740567  0.001213093  0.001757100  0.002764220  0.004157553;
0.000052080  0.001392067  0.001155927  0.001517920  0.001908893  0.002295547  0.003882547  0.004701120;
0.000356740  0.001475227  0.000969600  0.001911627  0.003466180  0.003516673  0.004207180  0.004477360;
0.001915233  0.001034687  0.001585680  0.001930367  0.002237273  0.002846040  0.004332327  0.003872000;
0.003726627  0.002722693  0.002572540  0.002523373  0.002708107  0.003281333  0.003654087  0.004934127;
0.015037607  0.011479927  0.009159580  0.008377720  0.008103813  0.008167760  0.009307100  0.008906027;
0.056371007  0.036370013  0.029310827  0.028868320  0.027531120  0.027664940  0.026556067  0.026637340;
0.229047847  0.149588360  0.123413240  0.110354107  0.102206900  0.100108547  0.098654247  0.106797240;

];

parallel_timings_with_cache_block_128 = [
    
0.000014933  0.000794860  0.000827713  0.000785073  0.001177453  0.001697780  0.003091493  0.003863873;
0.000038933  0.001103167  0.000584660  0.001621087  0.001246080  0.002287907  0.003598840  0.003923887;
0.000244067  0.000946727  0.000854300  0.001144980  0.001061927  0.002090500  0.002798853  0.003922200;
0.000968227  0.000820940  0.001096107  0.001061113  0.001269847  0.002652700  0.003538787  0.003801313;
0.003149233  0.002385220  0.002472407  0.002116260  0.002816380  0.003087807  0.003836580  0.005034987;
0.012219947  0.008817420  0.007611540  0.007464720  0.007060007  0.007518333  0.007935447  0.008797413;
0.060805267  0.040563073  0.031691933  0.028691707  0.027427967  0.027366567  0.028708773  0.029115820;
0.200006847  0.133841820  0.108478473  0.106037187  0.098513700  0.092656887  0.094730213  0.095034153;

];


parallel_timings_with_cache_block_256 = [
    
0 0 0 0 0 0 0 0;
0.000066247  0.001192013  0.000739293  0.001008553  0.001478493  0.002205253  0.003540647  0.003873067;
0.000228573  0.000943920  0.000836667  0.001146433  0.001252180  0.001562947  0.003318660  0.003613633;
0.001518547  0.002010180  0.001934667  0.002407813  0.001796820  0.003620260  0.003935980  0.003540513;
0.003927300  0.003172113  0.002822207  0.003310633  0.003561700  0.004063367  0.005498767  0.005726987;
0.013962973  0.010338733  0.009052520  0.007876207  0.008303033  0.007747527  0.009082227  0.009455540;
0.058288333  0.040827687  0.033835607  0.030575580  0.029909000  0.029483453  0.030098360  0.030032767;
0.191064733  0.123790313  0.109304340  0.102582947  0.099295653  0.101080647  0.096836287  0.095260767;

];


for i = 1:size(parallel_timings, 1)
    subplot(4,2,i); 
    hold on;
    grid on;
    grid minor;
%     xlim([1 16]);
%     ylim([0 6]);
    hold on;
    c_par_simd = plot(T, parallel_timings_simd(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);
    p_par_simd = plot(T, parallel_timings_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);
    
    c_par = plot(T, parallel_timings(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);
    p_par = plot(T, parallel_timings(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);
%     plot(T, parallel_timings_with_cache_block_64(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);
%     plot(T, parallel_timings_with_cache_block_128(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);
%     plot(T, parallel_timings_with_cache_block_256(i,:).^-1 * sequential_averages(i), 'LineWidth', 2);

%     legend('Speedup compared to sequential C implementation','Speedup compared to sequential Python implementation','Cache blocking (64)', 'Cache blocking (128)')

    hold off;
    title(sprintf('%d element Vector',matrix_sizes(i)),'interpreter','latex','fontsize',12);
    set(gca,'TickLabelInterpreter','latex')
end

hL = legend([c_par_simd, p_par_simd, c_par,p_par],{'Speedup with SIMD (relative to sequential C)','Speedup  with SIMD (relative to sequential Python)','Speedup without SIMD (relative to sequential C)','Speedup without SIMD (relative to sequential Python)'},'interpreter','latex','fontsize',12);
handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.XLabel.Interpreter = 'latex';


handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup','interpreter','latex','fontsize',16);
xlabel(handle,'Number of threads','interpreter','latex','fontsize',16);
title(handle,'\textbf{Speedup of parallelising Matrix-Vector multiplication with OpenMP threads \& SIMD support}','interpreter','latex','fontsize',12);
handle.Title.Position(2) = 1.05

