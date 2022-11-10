% clc; clear all;

  symulator_main1([90, 90], [100, 90],[100, 100] ,[90, 100],   1 , 10 , [90 , 90], [100 , 100], 0.5, 0); %wspolrzedne stacji referencyjnych, odleglosc miedzy punktami pomiarowymi, liczba pomiarow, wsp skrajnych rogów obszaru pomiarowego, wariancja, srednia
% symulator_main1([0, 0],  [30, 0], [20, 30], [10, 30], 1,10 , [1 , 1], [29 , 29], 0.2, 0);
function[] = symulator_main1(BS_1, BS_2, BS_3, BS_4, odleglosc, liczba_pomiarow, punkt_A, punkt_C, wariancja, srednia)

punkt_B = [punkt_C(1), punkt_A(2)];
punkt_D = [punkt_A(1), punkt_C(2)];
close all; %zamykanie figur

Chan_ori_v1_koniec_na_14a_s = (' Chan 14a');
Chan_ori_v1_koniec_na_14b_s = (' Chan 14b');
Chan_ori_v1_s = (' Chan'); 
Chan_ori_v2_koniec_na_14a_s = (' Chan 14a v2');
Chan_ori_v2_s = (' Chan v2');
Kkrrit_gradient_s = (' Modified TOA');

punkty_ref = [BS_1(1) BS_1(2); BS_2(1) BS_2(2); BS_3(1) BS_3(2); BS_4(1) BS_4(2)];

[wspolrzedne_punktow, liczba_punktow] = sprawdzenie(odleglosc, punkt_A, punkt_B, punkt_C, punkt_D);
[RMSE2d_Chan_ori_v1_koniec_na_14a, RMSE2d_Chan_ori_v1_koniec_na_14b, RMSE2d_Chan_ori_v1, RMSE2d_Chan_ori_v2_koniec_na_14a, RMSE2d_Chan_ori_v2, RMSE2d_Kkrrit_gradient,  A, B, C, D, E, F, G, H, I, J, K, L] = pomiary(BS_1, BS_2, BS_3, BS_4, wspolrzedne_punktow, liczba_punktow, liczba_pomiarow, wariancja, srednia);

%**********WYKRESY*********
kreslenie_dystrybuanty(RMSE2d_Chan_ori_v1_koniec_na_14a,RMSE2d_Chan_ori_v1_koniec_na_14b, RMSE2d_Chan_ori_v1, RMSE2d_Chan_ori_v2_koniec_na_14a, RMSE2d_Chan_ori_v2, RMSE2d_Kkrrit_gradient, liczba_punktow);
zlimits = kreslenie_RMSE3D(RMSE2d_Chan_ori_v1_koniec_na_14a, wspolrzedne_punktow, Chan_ori_v1_koniec_na_14a_s, [0,0] , punkty_ref);
zlimits = kreslenie_RMSE3D(RMSE2d_Chan_ori_v1_koniec_na_14b, wspolrzedne_punktow, Chan_ori_v1_koniec_na_14b_s, [0,0] , punkty_ref);
zlimits = kreslenie_RMSE3D(RMSE2d_Chan_ori_v1, wspolrzedne_punktow, Chan_ori_v1_s, [0,0], punkty_ref);
zlimits = kreslenie_RMSE3D(RMSE2d_Chan_ori_v2_koniec_na_14a, wspolrzedne_punktow, Chan_ori_v2_koniec_na_14a_s, [0,0], punkty_ref);
zlimits = kreslenie_RMSE3D(RMSE2d_Chan_ori_v2, wspolrzedne_punktow, Chan_ori_v2_s, [0,0], punkty_ref);
zlimits = kreslenie_RMSE3D(RMSE2d_Kkrrit_gradient, wspolrzedne_punktow, Kkrrit_gradient_s, [0,0], punkty_ref);


% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, A,B, punkt_A, punkt_B, punkt_C, punkt_D, Chan_ori_v1_koniec_na_14a_s);
% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, C,D, punkt_A, punkt_B, punkt_C, punkt_D, Chan_ori_v1_koniec_na_14b_s);
% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, E,F, punkt_A, punkt_B, punkt_C, punkt_D, Chan_ori_v1_s);
% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, G,H, punkt_A, punkt_B, punkt_C, punkt_D, Chan_ori_v2_koniec_na_14a_s);
% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, I,J, punkt_A, punkt_B, punkt_C, punkt_D, Chan_ori_v2_s);
% kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, K,L, punkt_A, punkt_B, punkt_C, punkt_D, Kkrrit_gradient_s);

end

function [RMSE2d_Chan_ori_v1_koniec_na_14a, RMSE2d_Chan_ori_v1_koniec_na_14b, RMSE2d_Chan_ori_v1, RMSE2d_Chan_ori_v2_koniec_na_14a, RMSE2d_Chan_ori_v2, RMSE2d_Kkrrit_gradient,  A, B, C, D, E, F, G, H, I, J, K, L] = pomiary(BS_1, BS_2, BS_3, BS_4, wspolrzedne_punktow, liczba_punktow, liczba_pomiarow, wariancja, srednia)

counterA = 0;
counterC = 0;
counterE = 0;
counterG = 0;
counterI = 0;
counterK = 0;

loss_Chan_14a = 0;
loss_Chan_14b = 0;
loss_Chan = 0;
loss_Chan_14a_v2 = 0;
lossChan_v2 = 0;
lossModified_TOA = 0;

inv_loss_Chan_14a = 0;
inv_loss_Chan = 0;
inv_loss_Chan_14a_v2 = 0;
inv_lossChan_v2 = 0;

counter_flag14av1 = 0;
counter_flagv1 = 0;
counter_flag14av2 = 0;
counter_flagv2 = 0;
for l = 1:liczba_punktow
            
          [pomiary_odleglosci] = wyznacz_zmierzone_odleglosci(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2),wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
    
   for i = 1:liczba_pomiarow
       
          [r] = gauss(wariancja, srednia);    %bï¿½ï¿½d gaussowski        
       
          [xChan_ori_v1_koniec_na_14a(i), yChan_ori_v1_koniec_na_14a(i), flag14av1] = Chan_ori_v1_koniec_na_14a(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r);    
          [xChan_ori_v1_koniec_na_14b(i), yChan_ori_v1_koniec_na_14b(i)] = Chan_ori_v1_koniec_na_14b(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r);
          [xChan_ori_v1(i), yChan_ori_v1(i), flagv1] = Chan_ori_v1(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r); 
          [xChan_ori_v2_koniec_na_14a(i), yChan_ori_v2_koniec_na_14a(i), flag14av2] = Chan_ori_v2_koniec_na_14a(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r);         
          [xChan_ori_v2(i), yChan_ori_v2(i),flagv2] = Chan_ori_v2(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r);   
          [xKkrrit_gradient(i), yKkrrit_gradient(i)] = Kkrrit_gradient(BS_1(1),BS_1(2),BS_2(1),BS_2(2),BS_3(1),BS_3(2),BS_4(1),BS_4(2), pomiary_odleglosci,r);
          
          counter_flag14av1 = flag14av1 + counter_flag14av1;
          counter_flagv1 = flagv1 + counter_flagv1;
          counter_flag14av2 = flag14av2 + counter_flag14av2;
          counter_flagv2 = flagv2 + counter_flagv2;
          
          mnoznik = 1;
          blad_gps = 15;
          
          if sqrt(((wspolrzedne_punktow(l,1) - xChan_ori_v1_koniec_na_14a(i))^2) + ((wspolrzedne_punktow(l,2) - yChan_ori_v1_koniec_na_14a(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)   
           A(i,l) = xChan_ori_v1_koniec_na_14a(i);
           B(i,l) = yChan_ori_v1_koniec_na_14a(i);
          else
              A(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              B(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterA = counterA + 1;
          end
          
          if sqrt(((wspolrzedne_punktow(l,1) - xChan_ori_v1_koniec_na_14b(i))^2) + ((wspolrzedne_punktow(l,2) - yChan_ori_v1_koniec_na_14b(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)
           C(i,l) = xChan_ori_v1_koniec_na_14b(i);
           D(i,l) = yChan_ori_v1_koniec_na_14b(i);
          else
              C(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              D(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterC = counterC + 1;
          end
          
          if sqrt(((wspolrzedne_punktow(l,1) - xChan_ori_v1(i))^2) + ((wspolrzedne_punktow(l,2) - yChan_ori_v1(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)
           E(i,l) = xChan_ori_v1(i);
           F(i,l) = yChan_ori_v1(i);
          else
              E(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              F(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterE = counterE + 1;
          end
          
          if sqrt(((wspolrzedne_punktow(l,1) - xChan_ori_v2_koniec_na_14a(i))^2) + ((wspolrzedne_punktow(l,2) - yChan_ori_v2_koniec_na_14a(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)
           G(i,l) = xChan_ori_v2_koniec_na_14a(i);
           H(i,l) = yChan_ori_v2_koniec_na_14a(i);
          else
              G(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              H(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterG = counterG + 1;
          end
          
          if sqrt(((wspolrzedne_punktow(l,1) - xChan_ori_v2(i))^2) + ((wspolrzedne_punktow(l,2) - yChan_ori_v2(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)
           I(i,l) = xChan_ori_v2(i);
           J(i,l) = yChan_ori_v2(i);
          else
              I(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              J(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterI = counterI + 1;
          end
          
          if sqrt(((wspolrzedne_punktow(l,1) - xKkrrit_gradient(i))^2) + ((wspolrzedne_punktow(l,2) - yKkrrit_gradient(i))^2)) <=  mnoznik* sqrt(((BS_1(1) - BS_2(1))^2) + ((BS_1(2) - BS_2(2)))^2)
           K(i,l) = xKkrrit_gradient(i);
           L(i,l) = yKkrrit_gradient(i);
          else
              K(i,l) = wspolrzedne_punktow(l,1) + blad_gps;
              L(i,l) = wspolrzedne_punktow(l,2) + blad_gps;
              counterK = counterK + 1;
          end
   end
      lenA = size(A,2);
      lenC = size(C,2);
      lenE = size(E,2);
      lenG = size(G,2);
      lenI = size(I,2);
      lenK = size(K,2);
      
      RMSE2d_Chan_ori_v1_koniec_na_14a(l) = RMSE(A(:,lenA),B(:,lenA),liczba_pomiarow,lenA, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
      RMSE2d_Chan_ori_v1_koniec_na_14b(l) = RMSE(C(:,lenC),D(:,lenC),liczba_pomiarow,lenC, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
      RMSE2d_Chan_ori_v1(l) = RMSE(E(:,lenE),F(:,lenE),liczba_pomiarow,lenE, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
      RMSE2d_Chan_ori_v2_koniec_na_14a(l) = RMSE(G(:,lenG),H(:,lenG),liczba_pomiarow,lenG, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
      RMSE2d_Chan_ori_v2(l) = RMSE(I(:,lenI),J(:,lenI),liczba_pomiarow,lenI, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
      RMSE2d_Kkrrit_gradient(l) = RMSE(K(:,lenK),L(:,lenK),liczba_pomiarow,lenK, wspolrzedne_punktow(l,1),wspolrzedne_punktow(l,2));
end

loss_Chan_14a = (counterA/(liczba_pomiarow*liczba_punktow))*100;
loss_Chan_14b = (counterC/(liczba_pomiarow*liczba_punktow))*100;
loss_Chan = (counterE/(liczba_pomiarow*liczba_punktow))*100;
loss_Chan_14a_v2 = (counterG/(liczba_pomiarow*liczba_punktow))*100;
lossChan_v2 = (counterI/(liczba_pomiarow*liczba_punktow))*100;
lossModified_TOA = (counterK/(liczba_pomiarow*liczba_punktow))*100;

inv_loss_Chan_14a = (counter_flag14av1/(liczba_pomiarow*liczba_punktow))*100;
inv_loss_Chan = (counter_flagv1/(liczba_pomiarow*liczba_punktow))*100;
inv_loss_Chan_14a_v2 = (counter_flag14av2/(liczba_pomiarow*liczba_punktow))*100;
inv_lossChan_v2 = (counter_flagv2/(liczba_pomiarow*liczba_punktow))*100;

names = {'loss_Chan_14a = ';'loss_Chan_14b ='; 'loss_Chan =';  'loss_Chan_14a_v2 = ';'lossChan_v2 =';'lossModified_TOA ='; 'inv_loss_Chan_14a =';'inv_loss_Chan = ';'inv_loss_Chan_14a_v2 =';'inv_lossChan_v2 ='};
values =[loss_Chan_14a; loss_Chan_14b; loss_Chan; loss_Chan_14a_v2; lossChan_v2; lossModified_TOA; inv_loss_Chan_14a ; inv_loss_Chan ; inv_loss_Chan_14a_v2; inv_lossChan_v2];
T = table(values, 'RowNames',names);

% T = table(['loss_Chan_14a = ', loss_Chan_14a; 'loss_Chan_14b =',loss_Chan_14b; 'loss_Chan =',loss_Chan ;    'loss_Chan_14a_v2 = ', loss_Chan_14a_v2;'lossChan_v2 =',lossChan_v2;'lossModified_TOA =',lossModified_TOA;    'inv_loss_Chan_14a =',inv_loss_Chan_14a ;'inv_loss_Chan = ',inv_loss_Chan ;'inv_loss_Chan_14a_v2 =',inv_loss_Chan_14a_v2;    'inv_lossChan_v2 =',inv_lossChan_v2]);
% Write data to text file
% writetable(T, 'Loss.txt','WriteRowNames',true);

end

function [wspolrzedne_punktow, liczba_punktow] = sprawdzenie(odleglosc,punkt_A, punkt_B, punkt_C, punkt_D)
spr_pkt_ob = 0;   
    while (spr_pkt_ob == 0)
        [wspolrzedne_punktow, liczba_punktow] = punkty_obszaru(odleglosc, punkt_A, punkt_B, punkt_C, punkt_D); %wygenerowanie punktï¿½w na obszarze
        [spr_pkt_ob] = sprawdz_punkty_obszaru (wspolrzedne_punktow, liczba_punktow);  %sprawdzenie czy punkty siï¿½ nakï¿½adajï¿½
    end
end

function [] = kreslenie_dystrybuanty(RMSE2d_Chan_ori_v1_koniec_na_14a, RMSE2d_Chan_ori_v1_koniec_na_14b, RMSE2d_Chan_ori_v1, RMSE2d_Chan_ori_v2_koniec_na_14a, RMSE2d_Chan_ori_v2, RMSE2d_Kkrrit_gradient,  liczba_punktow)
RMSE2d_Chan_ori_v1_koniec_na_14a_max = max(RMSE2d_Chan_ori_v1_koniec_na_14a);
RMSE2d_Chan_ori_v1_koniec_na_14b_max = max(RMSE2d_Chan_ori_v1_koniec_na_14b);
RMSE2d_Chan_ori_v1_max = max(RMSE2d_Chan_ori_v1);
RMSE2d_Chan_ori_v2_koniec_na_14a_max = max(RMSE2d_Chan_ori_v2_koniec_na_14a);
RMSE2d_Chan_ori_v2_max = max(RMSE2d_Chan_ori_v2);
RMSE2d_Kkrrit_gradient_max = max(RMSE2d_Kkrrit_gradient);

RMSE2d_matrix = [RMSE2d_Chan_ori_v1_koniec_na_14a_max, RMSE2d_Chan_ori_v1_koniec_na_14b_max, RMSE2d_Chan_ori_v1_max, RMSE2d_Chan_ori_v2_koniec_na_14a_max,RMSE2d_Chan_ori_v2_max, RMSE2d_Kkrrit_gradient_max];
             
RMSE2d_max = max(RMSE2d_matrix);

RMSE2d_max = round(RMSE2d_max + 0.5) + 0 ;


[x, dystrybuanta_Chan_ori_v1_koniec_na_14a] = dystrybuanta(RMSE2d_Chan_ori_v1_koniec_na_14a, liczba_punktow, RMSE2d_max);
[x, dystrybuanta_Chan_ori_v1_koniec_na_14b] = dystrybuanta(RMSE2d_Chan_ori_v1_koniec_na_14b, liczba_punktow, RMSE2d_max);
[x, dystrybuanta_Chan_ori_v1] = dystrybuanta(RMSE2d_Chan_ori_v1, liczba_punktow, RMSE2d_max);
[x, dystrybuanta_Chan_ori_v2_koniec_na_14a] = dystrybuanta(RMSE2d_Chan_ori_v2_koniec_na_14a, liczba_punktow, RMSE2d_max);
[x, dystrybuanta_Chan_ori_v2] = dystrybuanta(RMSE2d_Chan_ori_v2, liczba_punktow, RMSE2d_max);
[x, dystrybuanta_RMSE2d_Kkrrit_gradient] = dystrybuanta(RMSE2d_Kkrrit_gradient, liczba_punktow, RMSE2d_max);

figure();
plot(x,dystrybuanta_Chan_ori_v1_koniec_na_14a,'-or', 'MarkerSize',10);
hold on;
plot(x,dystrybuanta_Chan_ori_v1_koniec_na_14b,'-^g', 'MarkerSize',10);
plot(x,dystrybuanta_Chan_ori_v1,'-db', 'MarkerSize',10);
plot(x,dystrybuanta_Chan_ori_v2_koniec_na_14a,'-sm', 'MarkerSize',10);
plot(x,dystrybuanta_Chan_ori_v2,'-*c', 'MarkerSize',10);
plot(x,dystrybuanta_RMSE2d_Kkrrit_gradient,'-*b', 'MarkerSize',10);

title('Estymaty dystrybuanty RMSE')
xlabel('RMSE [m]') 
ylabel('Estymata dystrybuanty') 
legend ('Chan 14a v1', 'Chan 14b v1', 'Chan v1', 'Chan 14a v2', 'Chan v2', 'Modified TOA');
grid on;
hold off

end

function [zlimits] = kreslenie_RMSE3D(RMSE2d_Chan_ori_v1_koniec_na_14a, wspolrzedne_punktow, name, zlimits,  punkty_ref)

RMSE2d_TAC_transp = RMSE2d_Chan_ori_v1_koniec_na_14a';

tablica_do_wykresu_TAC(:,1) = RMSE2d_TAC_transp(:,1);
tablica_do_wykresu_TAC(:,2) = wspolrzedne_punktow(:,1);
tablica_do_wykresu_TAC(:,3) = wspolrzedne_punktow(:,2);

tablica_do_wykresu_psrt_TAC = sortrows(tablica_do_wykresu_TAC);

z = tablica_do_wykresu_psrt_TAC(:,1);
x = tablica_do_wykresu_psrt_TAC(:,2); 
y = tablica_do_wykresu_psrt_TAC(:,3);

SamplePerDim = 400;
X = linspace(min(x),max(x),SamplePerDim);
Y = linspace(min(y),max(y),SamplePerDim);
[X,Y] = ndgrid(X,Y);
F = scatteredInterpolant(x,y,z,'linear','none'); % scatteredInterpolant
Z = F(X,Y);
figure(),clf()
% scatter3(punkty_ref(:,1), y(:,1), z(:,1)+0.1, 'ko', 'linewidth',10)
scatter3 (punkty_ref(1,1), punkty_ref(1,2), 26, '^g', 'linewidth', 20);
hold on
scatter3 (punkty_ref(2,1), punkty_ref(2,2),26, '^r', 'linewidth', 20);
scatter3 (punkty_ref(3,1), punkty_ref(3,2),26, '^r', 'linewidth', 20);
scatter3 (punkty_ref(4,1), punkty_ref(4,2),26, '^r', 'linewidth', 20);
surf(X,Y,Z,'EdgeColor','none')
if(zlimits(1) == 0 && zlimits(2) == 0)
    zlimits = caxis;
else
    caxis(zlimits);
end
colormap(jet)
nazwa_colorbar = colorbar;
view(0,90)
title(nazwa_colorbar, 'RMSE [m]');
title(['RMSE 3D', name])
xlabel('X [m]') 
ylabel('Y [m]')
zlabel('RMSE [m]')


end

function [] = kreslenie_estymat(wspolrzedne_punktow, punkty_ref, liczba_punktow, X,Y, punkt_A, punkt_B, punkt_C, punkt_D, name )
figure();

hold on
plot (punkty_ref(:,1), punkty_ref(:,2), '^k', 'linewidth', 10);
    for i = 1:liczba_punktow
    plot (X(:,i),Y(:,i), '+', 'Color', [rand,rand,rand]);
    end
plot([punkt_A(1) punkt_B(1)],[punkt_A(2) punkt_B(2)],'r');
plot([punkt_B(1) punkt_C(1)],[punkt_B(2) punkt_C(2)],'r');
plot([punkt_C(1) punkt_D(1)],[punkt_C(2) punkt_D(2)],'r');
plot([punkt_D(1) punkt_A(1)],[punkt_D(2) punkt_A(2)],'r');

plot([punkty_ref(1,1),punkty_ref(2,1)], [punkty_ref(1,2),punkty_ref(2,2)],'b');
plot([punkty_ref(2,1),punkty_ref(3,1)], [punkty_ref(2,2),punkty_ref(3,2)],'b');
plot([punkty_ref(3,1),punkty_ref(4,1)], [punkty_ref(3,2),punkty_ref(4,2)],'b');
plot([punkty_ref(4,1),punkty_ref(1,1)], [punkty_ref(4,2),punkty_ref(1,2)],'b');

plot(wspolrzedne_punktow(:,1), wspolrzedne_punktow(:,2), 'sr','linewidth', 10, 'MarkerFaceColor', 'none');
title(['Estymaty po³o¿enia obiektów algorytmu',name])
xlabel('X [m]') 
ylabel('Y [m]') 

hold off
end

function [blad] = gauss(a,b)
rng('shuffle');
blad = a.*randn(1,4) + b;
end

function [x, ilosc_wartosci] = dystrybuanta(RMSE2d, liczba_punktow, max)
    
    for i = 1 : 0.1 : (1 + max)
    maksymalny_blad = i-0.9;
    k = round((i*10)-9);
    x(1)=0;
    
    if maksymalny_blad < 15
    x(k+1) = i-0.9;
    ilosc_wartosci(1)=0;
        ilosc_wartosci(k+1) = (sum(RMSE2d < (maksymalny_blad)))/ liczba_punktow;
    end
    end
    
    
   end

function [odleglosc] = wyznacz_zmierzone_odleglosci (x1_ref,y1_ref,x2_ref,y2_ref,x3_ref,y3_ref,x4_ref,y4_ref, wspolrzedne_punktow_x, wspolrzedne_punktow_y)

odleglosc(1) = sqrt(((wspolrzedne_punktow_x - x1_ref)^2) + ((wspolrzedne_punktow_y - y1_ref)^2)); 
odleglosc(2) = sqrt(((wspolrzedne_punktow_x - x2_ref)^2) + ((wspolrzedne_punktow_y - y2_ref)^2)); 
odleglosc(3) = sqrt(((wspolrzedne_punktow_x - x3_ref)^2) + ((wspolrzedne_punktow_y - y3_ref)^2));
odleglosc(4) = sqrt(((wspolrzedne_punktow_x - x4_ref)^2) + ((wspolrzedne_punktow_y - y4_ref)^2));

end

function [RMSE1] = RMSE(A,B,liczba_pomiarow, numer_punktu, wspolrzedne_punktow_x, wspolrzedne_punktow_y) 

RMSE1 = sqrt( sum((A-wspolrzedne_punktow_x).^2 + (B-wspolrzedne_punktow_y).^2)/liczba_pomiarow);

end

function [wspolrzedne_punktow, m] = punkty_obszaru(odleglosc,punkt_A, punkt_B, punkt_C, punkt_D)
    
dx = sqrt((punkt_A(1)-punkt_B(1))^2 + (punkt_A(2)-punkt_B(2))^2) ;
    dy = sqrt((punkt_A(1)-punkt_D(1))^2 + (punkt_D(2)-punkt_A(2))^2) ;
    krok_x = floor(dx / odleglosc)+1;
    krok_y = floor(dy / odleglosc)+1;
    wspolrzedne_punktow(1,1) = punkt_A(1);
    wspolrzedne_punktow(1,2) = punkt_A(2);
    
   m = krok_y*krok_x;
  
   w=0;
    for k=1:1:krok_y
        for i=1:1:krok_x  
            wspolrzedne_punktow(i+w,1) = punkt_A(1) + (i-1)*odleglosc;
            wspolrzedne_punktow(i+w,2) = punkt_A(2) + (k-1)*odleglosc;            
        end
        w=(numel( wspolrzedne_punktow)/2)+1;
         wspolrzedne_punktow(w,1) = punkt_A(1);
         wspolrzedne_punktow(w,2) = punkt_A(2)+k*odleglosc;
    
    end
    wspolrzedne_punktow= unique(wspolrzedne_punktow, 'rows', 'stable');
    wspolrzedne_punktow(end,:) = [];
end

function [spr_pkt_ob] = sprawdz_punkty_obszaru (wspolrzedne_punktow, liczba_punktow)

for x = 1:1:liczba_punktow 
    for y = 1:1:liczba_punktow 
        if (wspolrzedne_punktow(x,1) == wspolrzedne_punktow(y,1)) && (wspolrzedne_punktow(x,2) == wspolrzedne_punktow(y,2) && x ~= y)
            spr_pkt_ob = 0;
        else 
            spr_pkt_ob = 1;
        end
    end
end 
end

function [xChan, yChan, flag] = Chan_ori_v1(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

threshold = 10^(-10);
flag = 0;

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

x21 = x2 - x1;
x31 = x3 - x1;
x41 = x4 - x1;

y21 = y2 - y1;
y31 = y3 - y1;
y41 = y4 - y1;

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% ró¿nice odleg³oœci

r21 = r2 - r1;
r31 = r3 - r1;
r41 = r4 - r1;


%% wyliczenie r1

h = 0.5 * [r21^2 - K2 + K1;
           r31^2 - K3 + K1;
           r41^2 - K4 + K1 ];
       
Ga = - [x21, y21, r21;
        x31, y31, r31;
        x41, y41, r41 ];

Q = [1 0.5 0.5;
     0.5 1 0.5;
     0.5 0.5 1;];

za = ((Ga' * Q^-1 * Ga)^-1) * Ga' * Q^-1 * h;

x_approx_I = za(1);
y_approx_I = za(2);
r_approx_I = za(3);

%% iteracja nr 1

[x_iter_I, y_iter_I, r_iter_I, ~ , B, flag] = iteracja_v1(r21, r31, r41, Q, Ga, h, r_approx_I, za);

%% iteracja nr 2

[x_iter_II, y_iter_II, r_iter_II, psi, B, flag ] = iteracja_v1(r21, r31, r41, Q, Ga, h, r_iter_I, za);

%% wynik koñcowy
b_prim = [x_iter_II - x1, y_iter_II - y1, r_iter_II];
B_prim = diag(b_prim);

if(det(psi) < threshold)
    disp('psi jest nieodwracalna v2');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end

cov_za = (Ga'*psi^-1*Ga)^-1;
psi_prim = 4* B_prim * cov_za * B_prim;

if(det(psi_prim) < threshold)
    disp('psi_prim jest nieodwracalna v1');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end


Ga_prim = [1, 0;
           0, 1;
           1, 1;];
        
h_prim = [(x_iter_II-x1)^2;
          (y_iter_II-y1)^2;
          r_iter_II^2];

      
if(det(Ga_prim' * psi_prim^-1 * Ga_prim) < threshold)
    disp('Ga_prim * psi_prim^-1 * Ga_prim jest nieodwracalna v1');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end
      
      
za_prim = (Ga_prim' * psi_prim^-1 * Ga_prim)^-1 * Ga_prim' * psi_prim^-1 * h_prim;
za_prim_approx = (Ga_prim' * B_prim^-1 * Ga * Q^-1 * Ga * B_prim^-1 * Ga_prim) ^-1 * (Ga_prim' * B_prim^-1 * Ga * Q^-1 * Ga * B_prim^-1) * h_prim;

zp = sqrt(abs(za_prim)) + [x1; y1];
zp_approx = sqrt(za_prim_approx) + [x1; y1];

xChan = zp(1);
yChan = zp(2);

end

function [xChan, yChan] = Chan_ori_v1_koniec_na_14b(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

x21 = x2 - x1;
x31 = x3 - x1;
x41 = x4 - x1;

y21 = y2 - y1;
y31 = y3 - y1;
y41 = y4 - y1;

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% ró¿nice odleg³oœci

r21 = r2 - r1;
r31 = r3 - r1;
r41 = r4 - r1;


%% wyliczenie r1

h = 0.5 * [r21^2 - K2 + K1;
           r31^2 - K3 + K1;
           r41^2 - K4 + K1 ];
       
Ga = - [x21, y21, r21;
        x31, y31, r31;
        x41, y41, r41 ];

Q = [1 0.5 0.5;
     0.5 1 0.5;
     0.5 0.5 1;];
 
za = ((Ga' * Q^-1 * Ga)^-1) * Ga' * Q^-1 * h;

xChan = za(1);
yChan = za(2);
r_approx_I = za(3);

end

function [xChan, yChan, flag] = Chan_ori_v1_koniec_na_14a(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

x21 = x2 - x1;
x31 = x3 - x1;
x41 = x4 - x1;

y21 = y2 - y1;
y31 = y3 - y1;
y41 = y4 - y1;

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% ró¿nice odleg³oœci

r21 = r2 - r1;
r31 = r3 - r1;
r41 = r4 - r1;


%% wyliczenie r1

h = 0.5 * [r21^2 - K2 + K1;
           r31^2 - K3 + K1;
           r41^2 - K4 + K1 ];
       
Ga = - [x21, y21, r21;
        x31, y31, r31;
        x41, y41, r41 ];

Q = [1 0.5 0.5;
     0.5 1 0.5;
     0.5 0.5 1;];

za = ((Ga' * Q^-1 * Ga)^-1) * Ga' * Q^-1 * h;

x_approx_I = za(1);
y_approx_I = za(2);
r_approx_I = za(3);

%% iteracja nr 1

[x_iter_I, y_iter_I, r_iter_I, ~ , B, flag] = iteracja_v1(r21, r31, r41, Q, Ga, h, r_approx_I, za);

%% iteracja nr 2

[x_iter_II, y_iter_II, r_iter_II, psi, B, flag ] = iteracja_v1(r21, r31, r41, Q, Ga, h, r_iter_I, za);

xChan = x_iter_II;
yChan = y_iter_II;

end

function [xChan, yChan, flag] = Chan_ori_v2(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

threshold = 10^(-10);
flag = 0;

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

x21 = x2 - x1;
x31 = x3 - x1;
x41 = x4 - x1;

y21 = y2 - y1;
y31 = y3 - y1;
y41 = y4 - y1;

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% ró¿nice odleg³oœci

r21 = r2 - r1;
r31 = r3 - r1;
r41 = r4 - r1;


%% wyliczenie r1

h = 0.5 * [r21^2 - K2 + K1;
           r31^2 - K3 + K1;
           r41^2 - K4 + K1 ];
       
Ga = - [x21, y21, r21;
        x31, y31, r31;
        x41, y41, r41 ];

Q = [1 0.5 0.5;
     0.5 1 0.5;
     0.5 0.5 1;];

za = ((Ga' * Q^-1 * Ga)^-1) * Ga' * Q^-1 * h;

x_approx_I = za(1);
y_approx_I = za(2);
r_approx_I = za(3);

%% iteracja nr 1

[x_iter_I, y_iter_I, r_iter_I, ~ , B, flag] = iteracja_v2(r21, r31, r41, Q, Ga, h, x_approx_I, y_approx_I, x1, y1, za);

%% iteracja nr 2

[x_iter_II, y_iter_II, r_iter_II, psi, B, flag ] = iteracja_v2(r21, r31, r41, Q, Ga, h, x_iter_I, y_iter_I, x1, y1, za);

%% wynik koñcowy
b_prim = [x_iter_II - x1, y_iter_II - y1, r_iter_II];
B_prim = diag(b_prim);

if(det(psi) < threshold)
    disp('psi jest nieodwracalna v2');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end


cov_za = (Ga'*psi^-1*Ga)^-1;

psi_prim = 4* B_prim * cov_za * B_prim;


if(det(psi_prim) < threshold)
    disp('psi_prim jest nieodwracalna v2');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end


Ga_prim = [1, 0;
           0, 1;
           1, 1;];
        
h_prim = [(x_iter_II-x1)^2;
          (y_iter_II-y1)^2;
          r_iter_II^2];


if(det(Ga_prim' * psi_prim^-1 * Ga_prim) < threshold)
    disp('Ga_prim * psi_prim^-1 * Ga_prim jest nieodwracalna v2');
    flag = 1;
    xChan=za(1);
    yChan=za(2);
    return;
end
      
      
za_prim = (Ga_prim' * psi_prim^-1 * Ga_prim)^-1 * Ga_prim' * psi_prim^-1 * h_prim;
za_prim_approx = (Ga_prim' * B_prim^-1 * Ga * Q^-1 * Ga * B_prim^-1 * Ga_prim) ^-1 * (Ga_prim' * B_prim^-1 * Ga * Q^-1 * Ga * B_prim^-1) * h_prim;

zp = sqrt(abs(za_prim)) + [x1; y1];
% zp_approx = sqrt(za_prim_approx) + [x1; y1];
xChan = zp(1);
yChan = zp(2);

end

function [xChan, yChan, flag] = Chan_ori_v2_koniec_na_14a(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

x21 = x2 - x1;
x31 = x3 - x1;
x41 = x4 - x1;

y21 = y2 - y1;
y31 = y3 - y1;
y41 = y4 - y1;

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% ró¿nice odleg³oœci

r21 = r2 - r1;
r31 = r3 - r1;
r41 = r4 - r1;


%% wyliczenie r1

h = 0.5 * [r21^2 - K2 + K1;
           r31^2 - K3 + K1;
           r41^2 - K4 + K1 ];
       
Ga = - [x21, y21, r21;
        x31, y31, r31;
        x41, y41, r41 ];

Q = [1 0.5 0.5;
     0.5 1 0.5;
     0.5 0.5 1;];

za = ((Ga' * Q^-1 * Ga)^-1) * Ga' * Q^-1 * h;

x_approx_I = za(1);
y_approx_I = za(2);
r_approx_I = za(3);

%% iteracja nr 1

[x_iter_I, y_iter_I, r_iter_I, ~ , B, flag] = iteracja_v2(r21, r31, r41, Q, Ga, h, x_approx_I, y_approx_I, x1, y1, za);

%% iteracja nr 2

[x_iter_II, y_iter_II, r_iter_II, psi, B, flag ] = iteracja_v2(r21, r31, r41, Q, Ga, h, x_iter_I, y_iter_I, x1, y1, za);


xChan = x_iter_II;
yChan = y_iter_II;

end

function [x_iter, y_iter, r_iter, psi, B, flag] = iteracja_v1(r21, r31, r41, Q, Ga, h, r_iter, za)
flag = 0;
threshold = 10^(-10);

b = [r21 + r_iter, r31 + r_iter, r41 + r_iter];
B = diag(b);
    
psi = B * Q * B;

if(det(psi) < threshold)
    disp('psi jest nieodwracalna iter v1');
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end

if(det(Ga' * psi^-1 * Ga) < threshold)
    disp('Ga_prim * psi_prim^-1 * Ga_prim jest nieodwracalna iter v1');
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end
za = ((Ga' * psi^-1 * Ga)^-1) * Ga' * psi^-1 * h;
if(det(Ga' * psi^-1 * Ga) < threshold || det(psi)<threshold)
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end
if za(1) < 0
    za(1) = -za(1);
end

if za(2) < 0
    za(2) = -za(2);
end

x_iter = za(1);
y_iter = za(2);
r_iter = za(3);
end

function [x_iter, y_iter, r_iter, psi, B, flag] = iteracja_v2(r21, r31, r41, Q, Ga, h, x_approx, y_approx, x1, y1, za)

threshold = 10^(-10);
flag = 0;

r_iter = sqrt((x1 - x_approx )^2 + (y1 - y_approx)^2); %pomiary odl II wersja

b = [r21 + r_iter, r31 + r_iter, r41 + r_iter];
B = diag(b);
    
psi = B * Q * B;
if(det(psi) < threshold)
    disp('psi jest nieodwracalna iter v2');
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end

if(det(Ga' * psi^-1 * Ga) < threshold)
    disp('Ga_prim * psi_prim^-1 * Ga_prim jest nieodwracalna iter v2');
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end
za = ((Ga' * psi^-1 * Ga)^-1) * Ga' * psi^-1 * h;

if(det(Ga' * psi^-1 * Ga) < threshold || det(psi)<threshold)
    flag = 1;
    x_iter=za(1);
    y_iter=za(2);
    return;
end

if za(1) < 0
    za(1) = -za(1);
end

if za(2) < 0
    za(2) = -za(2);
end

x_iter = za(1);
y_iter = za(2);
r_iter = za(3);
end

function [xChan, yChan] = Kkrrit_gradient(x1, y1, x2, y2, x3, y3, x4, y4, pomiary_odleglosci,r)

r1 = pomiary_odleglosci(1) + r(1);
r2 = pomiary_odleglosci(2) + r(2);
r3 = pomiary_odleglosci(3) + r(3);
r4 = pomiary_odleglosci(4) + r(4);

K1 = x1^2 + y1^2;
K2 = x2^2 + y2^2;
K3 = x3^2 + y3^2;
K4 = x4^2 + y4^2;

%% 2013 v3

h = [r1^2 - K1;
     r2^2 - K2;
     r3^2 - K3;
     r4^2 - K4;];

Ga = [-2*x1, -2*y1, 1;
      -2*x2, -2*y2, 1;
      -2*x3, -2*y3, 1;
      -2*x4, -2*y4, 1];
 
q = [0.5, 0.5, 0.5, 0.5];
Q = diag(q);
   
za = (Ga' * Q^-1 * Ga)^-1 * (Ga' * Q^-1 * h);

  [x, y, za1] = iteracja_kkrrit ( x1, y1, x2, y2, x3, y3, x4, y4, r1, r2, r3, r4, za, Q, Ga, h);
% %% 2 iter
%  [x, y, za2] = iteracja_kkrrit ( x1, y1, x2, y2, x3, y3, x4, y4, r1, r2, r3, r4, za1, Q, Ga, h);
% %% 3 iter
%  [x, y, za3] = iteracja_kkrrit ( x1, y1, x2, y2, x3, y3, x4, r1, r2, r3, r4, y4, za2, Q, Ga, h);

% x = sqrt(abs(za(1)));
% y = sqrt(abs(za(2)));

xChan = x;
yChan = y;

end

function [x, y, za] = iteracja_kkrrit(x1, y1, x2, y2, x3, y3, x4, y4, r1, r2, r3, r4, za, Q, Ga, h )

B = [sqrt((x1 - za(1,1))^2 + (y1 - za(2,1))^2);
     sqrt((x2 - za(1,1))^2 + (y2 - za(2,1))^2);
     sqrt((x3 - za(1,1))^2 + (y3 - za(2,1))^2);
     sqrt((x4 - za(1,1))^2 + (y4 - za(2,1))^2)];
B = diag(B);
psi = 4* B * Q * B;

za = (Ga' * psi^-1 * Ga)^-1 * (Ga' * psi^-1 * h);

n = -[r1-B(1);
     r2-B(2);
     r3-B(3);
     r4-B(4);];

delta_za = (Ga' * psi^-1 * Ga)^-1 * Ga' * psi^-1 * B * n;
cov_za = (Ga' * psi^-1 * Ga)^-1;

b_prim = [za(1) - delta_za(1), za(2) - delta_za(2), 0.5];
B_prim = diag(b_prim);

psi_prim = 4 * B_prim * cov_za * B_prim;

Ga_prim = [1,0;
           0,1;
           1,1];
       
h_prim = [za(1,1)^2;
          za(2,1)^2;
          za(3,1)];
      
za = (Ga_prim' * psi_prim^-1 * Ga_prim)^-1 * (Ga_prim' * psi_prim^-1 * h_prim);

x = sqrt(abs(za(1)));
y = sqrt(abs(za(2)));

end