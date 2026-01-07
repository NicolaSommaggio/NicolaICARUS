import matplotlib.pyplot as plt
from collections import defaultdict
import random

# --- Scelte da tastiera ---
slice_da_disegnare = 1
vista_x = 2  # 1 = Induzione1, 2 = Induzione2, 3 = Collezione

# Colori fissi
colori_protoni = ["blue", "deepskyblue", "green", "orange", "darkgreen"]
indice_protone = 0
colori_assegnati = {}  # memorizza colore per nome-traccia base

# Funzione per colore random
def colore_random():
    return (random.random(), random.random(), random.random())

# Lettura file
with open("tracce_Run9435_1d_deconv.txt") as f:
    lines = f.readlines()

# Suddivisione in slice
slices = []
slice_corrente = []
for line in lines:
    if line.strip().startswith("#new slice"):
        slices.append(slice_corrente)
        slice_corrente = []
    else:
        slice_corrente.append(line.strip())
if slice_corrente:
    slices.append(slice_corrente)

# Controllo validità slice
if slice_da_disegnare > len(slices) or slice_da_disegnare < 1:
    raise ValueError("Slice non valida!")

dati = slices[slice_da_disegnare - 1]

# --- Parsing ---
tracce = defaultdict(list)
vertici = []

for line in dati:
    if not line or line.startswith("#"):
        continue

    if line.startswith("VERTEX"):
        vertici.append([float(x) for x in line.split()[1:]])
        continue

    campi = line.split()
    if len(campi) < 6:
        continue

    tempo = float(campi[0])
    induz1 = float(campi[1])
    induz2 = float(campi[2])
    coll = float(campi[3])
    nome_traccia = campi[4]
    tipo_particella = campi[5].lower()

    # Controllo se è una traccia sp
    is_sp = (len(campi) >= 7 and campi[6].lower() == "sp")

    #print(tempo,induz1,induz2,coll,nome_traccia,tipo_particella,is_sp)

    # Scelta vista X
    if vista_x == 1:
        x_coord = induz1
    elif vista_x == 2:
        x_coord = induz2
    elif vista_x == 3:
        x_coord = coll
    else:
        raise ValueError("Vista X non valida (scegli 1, 2 o 3)")

    # Salta punti con filo <= 0
    if x_coord <= 0:
        continue

    tracce[nome_traccia].append((tempo, x_coord, tipo_particella, is_sp))

# --- Plot ---
plt.figure(figsize=(10, 6))

for nome, punti in tracce.items():
    tempi = [p[0] for p in punti]
    fili = [p[1] for p in punti]
    tipo_particella = punti[0][2]
    is_sp = punti[0][3]

    #print(tempi[0],fili[0],tipo_particella,is_sp)
    """
    # Nome base (senza "sp")
    base_name = nome.replace("_sp", "").replace("sp", "")

    # Assegno colore al nome base
    if base_name not in colori_assegnati:
        if tipo_particella == "muon":
            colori_assegnati[base_name] = "red"
        elif tipo_particella == "proton":
            colori_assegnati[base_name] = colori_protoni[indice_protone % len(colori_protoni)]
            indice_protone += 1
        else:
            colori_assegnati[base_name] = colore_random()

    colore = colori_assegnati[base_name]
    """

    # Trasparenza per SP
    alpha_val = 0.7 if is_sp else 1.0

    # Etichetta sempre uguale alla base_name
    #etichetta = base_name
    etichetta = nome

    plt.plot(
        fili, tempi,
        linestyle="",
        marker="o",
        markersize=3,
        #color=colore,
        color=colore_random(),
        alpha=alpha_val,
        label=etichetta
    )

# --- Disegno vertici ---
for v in vertici:
    if vista_x == 1:
        vx = v[1]
    elif vista_x == 2:
        vx = v[2]
    elif vista_x == 3:
        vx = v[3]
    if vx > 0:
        plt.scatter(vx, v[0], color="black", s=90, marker="^",
                    facecolors='none', zorder=10, label="VERTEX")

# --- Cathode line ---
anode = 0
cathode_coord = 0
if v[4] > 0:
    anode = 359.33
    cathode_coord = 210
else:
    anode = -359.33
    cathode_coord = -210

cathode_time = abs(cathode_coord - anode) / 0.063 + 848
plt.axhline(y=cathode_time, color='yellow', linestyle='--', linewidth=1)

# --- Etichette assi ---
plt.xlabel(f"Wire number ({['Induction 1','Induction 2','Collection'][vista_x-1]})", fontsize=18)
plt.ylabel("Time", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

## --- Legenda pulita (una sola voce per traccia) ---
#handles, labels = plt.gca().get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys())
plt.legend()

# --- Salvataggio ---
plt.savefig('prova_event_display.pdf', format='pdf', bbox_inches='tight')
# plt.show()
