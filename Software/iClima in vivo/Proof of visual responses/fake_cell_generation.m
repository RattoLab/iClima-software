%%Script for the generation of moke cells (effect of spurious light on the
%%imaging)
G_VAL=mean(squeeze(mean(g(1:18,:,:),1,'omitnan')),1,'omitnan');
R_VAL=mean(squeeze(mean(r(1:18,:,:),1,'omitnan')),1,'omitnan');
mat_g=zeros(90,10,7);  %è una fake matrice di tutto il canale verde
%come se ogni cellula avesse un solo valore in tutti i trial e tutti i
%punti temporali. Nota: i valori che sono stati presi sono valori medi di
%ogni cellula del verde nei vari trial e nella baseline
mat_g(:,:,1)=G_VAL(1);
mat_g(:,:,2)=G_VAL(2);
mat_g(:,:,3)=G_VAL(3);
mat_g(:,:,4)=G_VAL(4);
mat_g(:,:,5)=G_VAL(5);
mat_g(:,:,6)=G_VAL(6);
mat_g(:,:,7)=G_VAL(7);

mat_r=zeros(90,10,7); %stessa cosa nel rosso.
mat_r(:,:,1)=R_VAL(1);
mat_r(:,:,2)=R_VAL(2);
mat_r(:,:,3)=R_VAL(3);
mat_r(:,:,4)=R_VAL(4);
mat_r(:,:,5)=R_VAL(5);
mat_r(:,:,6)=R_VAL(6);
mat_r(:,:,7)=R_VAL(7);

add_noise_g=mat_g+ch1_dark;
add_noise_r=mat_r+ch2_dark;

ratio_fake_cell=add_noise_g./add_noise_r;
for_origin=squeeze(mean(ratio_fake_cell,2,'omitnan'));

r0_fake=mean(for_origin(1:18,:),1,'omitnan');
fake_cells=-((for_origin-r0_fake)./r0_fake)*100;

mean(squeeze(mean(ch2_dark,2,'omitnan')),2)