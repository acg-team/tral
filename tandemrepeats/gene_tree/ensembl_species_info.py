''' translate ensembl_IDs, species (latin/standard) into each other'''

SPECIES_INFO = [
    ['Alpaca', 'Vicugna_pacos', 'ENSVPAP0', 'Eutheria', 8],
    ['Anole_lizard', 'Anolis_carolinensis', 'ENSACAP0', 'Amniota', 11],
    ['Armadillo', 'Dasypus_novemcinctus', 'ENSDNOP0', 'Eutheria', 8],
    ['Bushbaby', 'Otolemur_garnettii', 'ENSOGAP0', 'Primates', 6],
    ['Ciona_intestinalis', 'Ciona_intestinalis', 'ENSCINP0', 'Chordata', 16],
    ['Ciona_savignyi', 'Ciona_savignyi', 'ENSCSAVP0', 'Chordata', 16],
    ['Cat', 'Felis_catus', 'ENSFCAP0', 'Eutheria', 8],
    ['Chicken', 'Gallus_gallus', 'ENSGALP0', 'Amniota', 11],
    ['Chimpanzee', 'Pan_troglodytes', 'ENSPTRP0', 'Homininae', 0],
    ['Chinese_softshell_turtle', 'Pelodiscus_sinensis', 'ENSPSIP0', 'Amniota', 11],
    ['Cod', 'Gadus_morhua', 'ENSGMOP0', 'Euteleostomi', 14],
    ['Coelacanth', 'Latimeria_chalumnae', 'ENSLACP0', 'Sarcopterygii', 13],
    ['Cow', 'Bos_taurus', 'ENSBTAP0', 'Eutheria', 8],
    ['Dog', 'Canis_lupus_familiaris', 'ENSCAFP0', 'Eutheria', 8],
    ['Dolphin', 'Tursiops_truncatus', 'ENSTTRP0', 'Eutheria', 8],
    ['Elephant', 'Loxodonta_africana', 'ENSLAFP0', 'Eutheria', 8],
    ['Ferret', 'Mustela_putorius_furo', 'ENSMPUP0', 'Eutheria', 8],
    ['Fruitfly', 'Drosophila_melanogaster', 'FBpp', 'Bilateria', 17],
    ['Fugu', 'Takifugu_rubripes', 'ENSTRUP0', 'Euteleostomi', 14],
    ['Gibbon', 'Nomascus_leucogenys', 'ENSNLEP0', 'Hominoidea', 2],
    ['Gorilla', 'Gorilla_gorilla_gorilla', 'ENSGGOP0', 'Homininae', 0],
    ['Guinea_Pig', 'Cavia_porcellus', 'ENSCPOP0', 'Euarchontoglires', 7],
    ['Hedgehog', 'Erinaceus_europaeus', 'ENSEEUP0', 'Eutheria', 8],
    ['Horse', 'Equus_caballus', 'ENSECAP0', 'Eutheria', 8],
    ['Human', 'Homo_sapiens', 'ENSP0', 'Homo_sapiens', -1],
    ['Hyrax', 'Procavia_capensis', 'ENSPCAP0', 'Eutheria', 8],
    ['Kangaroo_rat', 'Dipodomys_ordii', 'ENSDORP0', 'Euarchontoglires', 7],
    ['Lamprey', 'Petromyzon_marinus', 'ENSPMAP0', 'Vertebrata', 15],
    ['Lesser_hedgehog_tenrec', 'Echinops_telfairi', 'ENSETEP0', 'Eutheria', 8],
    ['Macaque', 'Macaca_mulatta', 'ENSMMUP0', 'Catarrhini', 3],
    ['Marmoset', 'Callithrix_jacchus', 'ENSCJAP0', 'Simiiformes', 4],
    ['Medaka', 'Oryzias_latipes', 'ENSORLP0', 'Euteleostomi', 14],
    ['Megabat', 'Pteropus_vampyrus', 'ENSPVAP0', 'Eutheria', 8],
    ['Microbat', 'Myotis_lucifugus', 'ENSMLUP0', 'Eutheria', 8],
    ['Mouse', 'Mus_musculus', 'ENSMUSP0', 'Euarchontoglires', 7],
    ['Mouse_Lemur', 'Microcebus_murinus', 'ENSMICP0', 'Primates', 6],
    ['Opossum', 'Monodelphis_domestica', 'ENSMODP0', 'Theria', 9],
    ['Orangutan', 'Pongo_abelii', 'ENSPPYP0', 'Hominidae', 1],
    ['Panda', 'Ailuropoda_melanoleuca', 'ENSAMEP0', 'Eutheria', 8],
    ['Pig', 'Sus_scrofa', 'ENSSSCP0', 'Eutheria', 8],
    ['Pika', 'Ochotona_princeps', 'ENSOPRP0', 'Euarchontoglires', 7],
    ['Platyfish', 'Xiphophorus_maculatus', 'ENSXMAP0', 'Euteleostomi', 14],
    ['Platypus', 'Ornithorhynchus_anatinus', 'ENSOANP0', 'Mammalia', 10],
    ['Rabbit', 'Oryctolagus_cuniculus', 'ENSOCUP0', 'Euarchontoglires', 7],
    ['Rat', 'Rattus_norvegicus', 'ENSRNOP0', 'Euarchontoglires', 7],
    ['Roundworm','Caenorhabditis_elegans','XXXXXXX', 'Bilateria', 17],
    ['Shrew', 'Sorex_araneus', 'ENSSARP0', 'Eutheria', 8],
    ['Sloth', 'Choloepus_hoffmanni', 'ENSCHOP0', 'Eutheria', 8],
    ['Squirrel', 'Ictidomys_tridecemlineatus', 'ENSSTOP0', 'Euarchontoglires', 7],
    ['Stickleback', 'Gasterosteus_aculeatus', 'ENSGACP0', 'Euteleostomi', 14],
    ['Tarsier', 'Tarsius_syrichta', 'ENSTSYP0', 'Haplorrhini', 5],
    ['Tasmanian_devil', 'Sarcophilus_harrisii', 'ENSSHAP0', 'Theria', 9],
    ['Tetraodon', 'Tetraodon_nigroviridis', 'ENSTNIP0', 'Euteleostomi', 14],
    ['Tilapia', 'Oreochromis_niloticus', 'ENSONIP0', 'Euteleostomi', 14],
    ['Tree_Shrew', 'Tupaia_belangeri', 'ENSTBEP0', 'Euarchontoglires', 7],
    ['Turkey', 'Meleagris_gallopavo', 'ENSMGAP0', 'Amniota', 11],
    ['Wallaby', 'Macropus_eugenii', 'ENSMEUP0', 'Theria', 9],
    ['Xenopus', 'Xenopus_tropicalis', 'ENSXETP0', 'Tetrapoda', 12],
    ['Zebra_Finch', 'Taeniopygia_guttata', 'ENSTGUP0', 'Amniota', 11],
    ['Zebrafish', 'Danio_rerio', 'ENSDARP0', 'Euteleostomi', 14],
    ["Yeast",'Saccharomyces_cerevisiae','XXXXXXX', 'Opisthokonta', 18],
]

order = {'english':0, 'official':1, 'ensembl_ID':2, 'clade_name':3, 'clade_distance':4}

# ["Yeast","Saccharomyces_cerevisiae","???"]
YEAST = ['YI', 'Q0', 'YK', 'YJ', 'YM', 'YL', 'YO', 'YH', 'YA', 'YC', 'YB', 'YE', 'YD', 'YG', 'YF', 'YN', 'YP']
# ["Caenorhabditis_elegans","Caenorhabditis_elegans","???"]
WORM = ['2L', '2R', '3R', '4R', '6R', 'AC', 'AH', 'B0', 'BE', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'CC', 'CD', 'CE', 'cT', 'D1', 'D2', 'DC', 'DH', 'DY', 'E0', 'EE', 'EG', 'E_', 'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'H0', 'H1', 'H2', 'H3', 'H4', 'JC', 'K0', 'K1', 'LL', 'M0', 'M1', 'M2', 'M4', 'M5', 'M6', 'M7', 'M8', 'MT', 'PA', 'PD', 'R0', 'R1', 'R3', 'R5', 'R7', 'R9', 'SS', 'T0', 'T1', 'T2', 'VB', 'VC', 'VF', 'VH', 'VK', 'VM', 'VW', 'VY', 'VZ', 'W0', 'W1', 'Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9', 'ZC', 'ZK']

CLADE_NAME_TO_DISTANCE = {"Opisthokonta":18, "Bilateria":17}

def translate(translate_from, translate_to, text):

    if not (translate_from in order and translate_to in order):
        print('translate_from: {0} or translate_to: {1} not in order'.format(translate_from,translate_to))
        raise ValueError('translate_from: {0} or translate_to: {1} not in order'.format(translate_from,translate_to))
        return None
    
    if translate_to == 'ensembl_ID':
        if text in ['Yeast', 'Saccharomyces_cerevisiae']:
            return YEAST
        elif text in ['Yeast', 'Saccharomyces_cerevisiae']:
            return WORM
    elif translate_from == 'ensembl_ID':
        for iY in YEAST:
            if text.startswith(iY):
                if translate_to == 'english':
                    return "Yeast"
                elif translate_to == 'clade_distance':
                    return CLADE_NAME_TO_DISTANCE["Opisthokonta"]
                elif translate_to == 'clade_name':
                    return "Opisthokonta"
                else:
                    return "Saccharomyces_cerevisiae"
        for iW in WORM:
            if text.startswith(iW):
                if translate_to == 'clade_distance':
                    return CLADE_NAME_TO_DISTANCE["Bilateria"]
                elif translate_to == 'clade_name':
                    return "Bilateria"
                elif translate_to == 'english':
                    return "Roundworm"
                else:
                    return "Caenorhabditis_elegans"
    for species in SPECIES_INFO:
        if translate_from == 'ensembl_ID':
            if text.startswith(species[order[translate_from]]):
                return species[order[translate_to]]
        else:
            if text == species[order[translate_from]]:
                return species[order[translate_to]]
    else:
        print('text {0} not in translate_from: {1}'.format(text,translate_from))
        raise ValueError('text {0} not in translate_from: {1}'.format(text,translate_from))
        return None
