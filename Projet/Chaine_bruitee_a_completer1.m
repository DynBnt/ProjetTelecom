%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETUDE MODULATEURDEMODULATEUR SUR CANAL AWGN
% BENOIT Doryan, Mai 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES GENERAUX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = [1 0 0 0 1 0 1;
     0 1 0 0 1 1 1;
     0 0 1 0 1 1 0;
     0 0 0 1 0 1 1];

% Fonction de codage de Hamming
function c = codage_hamming(bits_reshape,G)
    % Matrice génératrice pour un code de Hamming (7,4)
    
    c = reshape(mod(bits_reshape * G,2),1,[]);
end

function bits_decodes = decode_bits_dur(bits_codes,G)
% decode_bits  Décodage des bits préalablement codés
% 
%   bits_codes : vecteur ligne 1×2100 de bits codés
%   G          : matrice de génération (4×7)
%
%   bits_decodes : vecteur ligne 1×1200 des bits d'information décodés

    % Regrouper en mots codés de 7 bits : M est 300×7
    M = reshape(bits_codes, 7, []).';  % reshape puis transpose
    
    % Construire la matrice b (16×4)
    b = de2bi(0:15, 4, 'left-msb');  % 16×4

    % Calculer le dictionnaire D = b * G mod 2  (16×7)
    D = mod(b * G, 2);
    
    %On crée la matrice U qui va accueillir les nouveaux bits décodés
    U = zeros(size(M,1), 4);

    % Pour chaque mot codé M_i, trouver le mot d'info le plus proche
    for i = 1:size(M,1)
        % Réplication de M_i en 16×7
        Mi_rep = repmat(M(i,:), 16, 1);

        % XOR avec le dictionnaire et matrice de M_i
        X = xor(D, Mi_rep);

        % Distance : somme des 1 ligne par ligne : matrice 16×1
        distances = sum(X, 2);

        % Indice du dictionnaire minimal
        [~, j] = min(distances);

        % Récupérer les 4 bits d'info correspondants
        U(i, :) = b(j, :);
    end

    % Mettre U en vecteur ligne 1×1200
    bits_decodes = reshape(U.', 1, []);  % transpose puis reshape

end


Fe=24000;       % Fréquence d'échantillonnage
Te=1/Fe;        % Période d'échantillonnage
Rb=3000;        % Débit binaire souhaité
N=1200;         % Nombre de bits générés

M=2;            % Ordre de la modulation (binaire à moyenne nulle)
Rs=Rb/log2(M);  % Débit symbole
Ns=Fe/Rs;       % Facteur de suréchantillonnage

% Tableau des valeurs de SNR par bit souhaité à l'entrée du récepteur en dB
tab_Eb_N0_dB=9:15; 
% Passage au SNR en linéaire
tab_Eb_N0=10.^(tab_Eb_N0_dB/10);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUCLE SUR LES NIVEAUX DE EbN0 A TESTER POUR OBTENTION DU TES ET DU TEB
% SIMULES DE LA CHAINE IMPLANTEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indice_bruit=1:length(tab_Eb_N0_dB)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VALEUR DE EbN0 TESTEE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eb_N0_dB=tab_Eb_N0_dB(indice_bruit);
    Eb_N0=tab_Eb_N0(indice_bruit);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALISATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nb_erreurs=0;   % Variable permettant de compter le nombre d'erreurs cumulées
    nb_cumul=0;     % Variables permettant de compter le nombre de cumuls réalisés
    TES=0;          % Initialisation du taux d'erreur symbole pour le cumul
    TEB=0;          % Initialisation du taux d'erreur binaire pour le cumul

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BOUCLE POUR PRECISION TEB MESURE  COMPTAGE NOMBRE ERREURS
    % (voir annexe texte TP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(nb_erreurs<1000)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATION DE L'INFORMATION BINAIRE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits=randi([0,1],1,N);
        bits_reshape = reshape(bits,N/4,4);
        bits_codees = codage_hamming(bits_reshape,G); %Vecteur des bits ( 1,2100 )
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        symboles = 2*bits_codees - 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SURECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        somme_Diracs_ponderes=kron(symboles,[1 zeros(1,Ns-1)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILTRAGE DE MISE EN FORME 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Génération de la réponse impulsionnelle du filtre de mise en forme
        h = ones(1,Ns);
        % Filtrage de mise en forme
        Signal_emis=filter(h,1,somme_Diracs_ponderes);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CANAL DE PROPAGATION AWGN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calcul de la puissance du signal émis
        P_signal = mean(Signal_emis.^2);
        % Calcul de la puissance du bruit à ajouter pour obtenir EbN0 souhaité
        P_bruit = P_signal  * Ns / (2*log2(M)*Eb_N0); 
        % Génération du bruit gaussien
        Bruit=sqrt(P_bruit)*randn(1,length(Signal_emis)); 
        % Ajout du bruit canal au signal émis
        Signal_recu=Signal_emis+Bruit;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILTRAGE DE RECEPTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hr = ones(1,Ns);
        Signal_recu_filtre=filter(hr,1,Signal_recu);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choix de n0
        n0 = Ns;
        % Echantillonnage à n0+mNs
        Signal_echantillonne=Signal_recu_filtre(n0:Ns:end);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DECISIONS SUR LES SYMBOLES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        symboles_recus = 2*(Signal_echantillonne>0) - 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCUL DU TAUX D'ERREUR SYMBOLE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TES = TES + sum(symboles_recus ~= symboles) / length(symboles);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEMAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits_recus = (symboles_recus + 1)/2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DECODAGE DUR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits_decodes = decode_bits_dur(bits_recus,G);
        bits_recus = bits_decodes;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CALCUL DU TAUX D'ERREUR BINAIRE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TEB = TEB + sum(bits_recus ~= bits) / length(bits);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CUMUL DU NOMBRE D'ERREURS ET NOMBRE DE CUMULS REALISES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb_erreurs = nb_erreurs + sum(bits_recus ~= bits);
        nb_cumul = nb_cumul + 1;

    end  % fin boucle sur comptage nombre d'erreurs

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DU TAUX D'ERREUR SYMBOLE ET DU TAUX D'ERREUR BINAIRE POUR LA
    % VALEUR TESTEE DE EbN0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TES_simule(indice_bruit)=TES/nb_cumul;       
    TEB_simule(indice_bruit)=TEB/nb_cumul; 

     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIAGRAMME DE L'OEIL EN SORTIE DU FILTRE DE RECEPTION AVEC BRUIT
    % TRACE POUR CHAQUE VALEUR DE EbN0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oeil=reshape(Signal_recu_filtre,Ns,length(Signal_recu_filtre)/Ns);
    figure
    plot(oeil)
    title(['Tracé du diagramme de loeil en sortie du filtre de réception pour E_bN_0 = ' num2str(Eb_N0_dB) ' dB'])

end  % fin boucle sur les valeurs testées de EbN0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DU TES ET DU TEB THEORIQUE DE LA CHAINE IMPLANTEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = sqrt(6*log2(M)/(M^2-1) * tab_Eb_N0);
TES_THEO = 2*(M-1)/M * qfunc(alpha);
TEB_THEO = 2*(M-1)/(M*log2(M)) * qfunc(alpha);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACES DES TES ET TEB OBTENUS EN FONCTION DE EbN0
% COMPARAISON AVEC LES TES et TEBs THEORIQUES DE LA CHAINE IMPLANTEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(tab_Eb_N0_dB, TES_THEO,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TES_simule,'b-o')
legend('TES théorique','TES simulé')
xlabel('E_bN_0 (dB)')
ylabel('TES')
title('TES M-PAM sur canal AWGN')

figure
semilogy(tab_Eb_N0_dB, TEB_THEO,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TEB_simule,'b-o')
legend('TEB théorique','TEB simulé')
xlabel('E_bN_0 (dB)')
ylabel('TEB')
title('TEB M-PAM sur canal AWGN')
