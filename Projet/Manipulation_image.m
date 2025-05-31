%Lecture de l'image
image = imread('dcode-image.png');
%Visualisation
figure
imshow(image)
%Transformation de l'image en un train binaire
vect_image=reshape(image,1,size(image,1)*size(image,2));
mat_image_binaire=de2bi(vect_image);
bits=double(reshape(mat_image_binaire,1,size(mat_image_binaire,1)*size(mat_image_binaire,2)));

%Reconstruction de l'image à partir de la suite binaire
mat_image_binaire_retrouvee=reshape(bits,211*300,8);
mat_image_decimal_retrouvee=bi2de(mat_image_binaire_retrouvee);
image_retrouvee=reshape(mat_image_decimal_retrouvee,211,300);
%Visualisation
figure
imshow(uint8(image_retrouvee))
