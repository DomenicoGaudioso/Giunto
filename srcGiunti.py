import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
#import pyvista as pv
#from pyvista.plotting import Plotter

# Titolo dell'app
st.title("Verifica giunzione metallica")

# Input dei parametri della giunzione
propB = {"M12": {"d": 12, "Ares": 84.3, "k":8, "e":21.9},
         "M14": {"d": 14, "Ares": 115, "k":9, "e":25.4},
         "M16": {"d": 16, "Ares": 157, "k":10, "e":27.7},
         "M18": {"d": 18, "Ares": 192, "k":11, "e":30},
         "M20": {"d": 20, "Ares": 245, "k":13, "e":34.6},
         "M22": {"d": 22, "Ares": 303, "k":14, "e":37},
         "M24": {"d": 24, "Ares": 353, "k":15, "e":41.6},
         "M27": {"d": 27, "Ares": 459, "k":17, "e":48},
         "M30": {"d": 30, "Ares": 561, "k":19, "e":53.1},
         "M36": {"d": 36, "Ares": 817, "k":23, "e":63.5},
         }

classeB = {"4.6": {"fyb": 240, "fub": 400},
         "4.8": {"fyb": 320, "fub": 400},
         "5.6": {"fyb": 300, "fub": 500},
         "5.8": {"fyb": 400, "fub": 500},
         "6.8": {"fyb": 480, "fub": 600},
         "8.8": {"fyb": 640, "fub": 800},
         "10.9": {"fyb": 900, "fub": 1000},
         }

acciaio = { "S235": {"fyk": 235, "ftk": 360,"fyk_40": 215, "ftk_40": 360},
            "S275": {"fyk": 275, "ftk": 430,"fyk_40": 255, "ftk_40": 410},
            "S355": {"fyk": 355, "ftk": 510,"fyk_40": 335, "ftk_40": 470},
            "S450": {"fyk": 440, "ftk": 550,"fyk_40": 420, "ftk_40": 550},
            
            "S275 N/NL": {"fyk": 275, "ftk": 390,"fyk_40": 255, "ftk_40": 370},
            "S355 N/NL": {"fyk": 355, "ftk": 490,"fyk_40": 335, "ftk_40": 470},
            "S420 N/NL": {"fyk": 420, "ftk": 520,"fyk_40": 390, "ftk_40": 520},
            "S460 N/NL": {"fyk": 460, "ftk": 540,"fyk_40": 430, "ftk_40": 540},
            
            "S275 M/ML": {"fyk": 275, "ftk": 370,"fyk_40": 255, "ftk_40": 360},
            "S355 M/ML": {"fyk": 355, "ftk": 470,"fyk_40": 335, "ftk_40": 450},
            "S420 M/ML": {"fyk": 420, "ftk": 520,"fyk_40": 390, "ftk_40": 500},
            "S460 M/ML": {"fyk": 460, "ftk": 540,"fyk_40": 430, "ftk_40": 530},
            "S460 Q/QL/QL1": {"fyk": 460, "ftk": 570,"fyk_40": 440, "ftk_40": 580},
            
            "S235 W": {"fyk": 235, "ftk": 360,"fyk_40": 215, "ftk_40": 340},
            "S355 W": {"fyk": 355, "ftk": 510,"fyk_40": 335, "ftk_40": 490},
         }


mu = {"superifici sabbiate neccanicamente o a graniglia, esenti da incrostazioni di ruggine e da vaiolature": 0.5,
      "superifici sabbiate neccanicamente o a graniglia, e verniciate a spruzzo con prodotti a base di alluminio o di zinco": 0.4,
      "superifici pulite mediante spazzolatura o alla fiamma, esenti da incrostazioni di ruggine": 0.3,
      "superifici non trattate": 0.2,
      }

gm0 = 1.05
gm1 = 1.1
gm2 = 1.25
gm3_slu = 1.25
gm3_sle = 1.1
gm2 = 1.25
gm6 = 1.0
gm7_1 = 1.0
gm7_2 = 1.10

# Selezione delle chiavi
st.sidebar.header("Parametri della Giunzione")

selected_propB_key = st.sidebar.selectbox("Seleziona il tipo di bullone (propB):", options=list(propB.keys()))
selected_classeB_key = st.sidebar.selectbox("Seleziona la classe del bullone (classeB):", options=list(classeB.keys()))
selected_acciaio_key = st.sidebar.selectbox("Seleziona acciao per la piastra:", options=list(acciaio.keys()))

tp = st.sidebar.number_input(key="tp",label= "spessore piastra (mm):", min_value=1, value=10, step=1)
npt = st.sidebar.number_input(key="pt",label= "piani di taglio:", min_value=1, value=1, step=1)

e1 = st.sidebar.number_input(key="e1", label="e1: distanza del bullone esterno dal bordo in direzione della forza (mm):", min_value=1, value=30, step=10)
p1 = st.sidebar.number_input(key="p1", label="p1: distanza tra bulloni interni in direzione della forza (mm):", min_value=0, value=60, step=10)
e2 = st.sidebar.number_input(key="e2", label="e2: distanza del bullone esterno dal bordo in direzione ortogonale alla forza (mm):", min_value=1, value=30, step=10)
p2 = st.sidebar.number_input(key="p2", label="p2: distanza tra bulloni interni in direzione ortogonali alla forza (mm):", min_value=0, value=60, step=10)
nb = st.sidebar.number_input(key="nb", label="nb: numero bulloni:", min_value=1, value=10, step=1)
nf = st.sidebar.number_input(key="nf", label="nf: numero file:", min_value=1, value=2, step=1)

db = selected_propB_value = propB[selected_propB_key]["d"]


# Sollecitazioni
V_slu = st.sidebar.number_input("Taglio SLU(kN):", min_value=0.0, value=150.0, step=1.0)
T_slu = st.sidebar.number_input("Trazione SLU(kN):", min_value=0.0, value=10.0, step=1.0)
#M_slu = st.sidebar.number_input("Momento flettente SLU(kNm):", min_value=0.0, value=50.0, step=1.0)

V_sle = st.sidebar.number_input("Taglio SLE (kN):", min_value=0.0, value=80.0, step=1.0)
T_sle = st.sidebar.number_input("Trazione SLE(kN):", min_value=0.0, value=0.0, step=1.0)
#M_sle = st.sidebar.number_input("Momento flettente SLE (kNm):", min_value=0.0, value=50.0, step=1.0)

nx = int(nb/nf) #numero di bulloni per fila
lpx = e1*2+p1*(nx-1)
lpy = e2*2+p2*(nf-1)


# Recupero dei valori selezionati
selected_propB_value = propB[selected_propB_key]
selected_classeB_value = classeB[selected_classeB_key]
selected_acciaio_value = acciaio[selected_acciaio_key]


d0 = db+1 if db<=20 else db+1.5 #diametro foro
#distanze minime
e1_min = 1.2*d0
e2_min = 1.2*d0
p1_min = 2.2*d0
p2_min = 2.4*d0
#distanze massime
tp_min = tp
e1_max = 4*tp_min + 40 
e2_max = 4*tp_min + 40 
p1_max = min(14*tp_min, 200)
p2_max = min(14*tp_min, 200)

# Disegno i bulloni
Bulloni = {} #coordinate dei bulloni
b = 0
for i in range(0,nf):
    for j in range(0,nx):
        b = b+1
        x = e1 + j*p1
        y = e2 + i*p2
        Bulloni[b] = [x, y]

#---------------------------------------------------------------------------#
# Funzione per creare un cilindro
def create_cylinder(x_center, y_center, z_bottom, height, radius, resolution=20):
    theta = np.linspace(0, 2 * np.pi, resolution)
    z = z_bottom  # Altezza del cilindro
    theta_grid, z_grid = np.meshgrid(theta, z)
    x = radius * np.cos(theta_grid) + x_center
    y = radius * np.sin(theta_grid) + y_center
    return x.flatten(), y.flatten(), z_grid.flatten()

# Creazione del parallelepipedo (piastra)
vertices = [
    [0, 0, 0], [lpx, 0, 0], [lpx, lpy, 0], [0, lpy, 0],  # Base inferiore
    [0, 0, tp], [lpx, 0, tp], [lpx, lpy, tp], [0, lpy, tp],  # Base superiore
]
faces = [
    [0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
    [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7],
]
x, y, z = zip(*vertices)
i, j, k = [], [], []

for face in faces:
    i.extend([face[0], face[0], face[1], face[1], face[2], face[2]])
    j.extend([face[1], face[3], face[2], face[0], face[3], face[0]])
    k.extend([face[3], face[2], face[3], face[2], face[0], face[1]])

fig = go.Figure()

# Aggiunta della piastra
fig.add_trace(
    go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k,
        color='lightblue', opacity=0.5, name="Piastra"
    )
)

dd = propB[selected_propB_key]["e"]
hd = propB[selected_propB_key]["k"]

# Aggiunta dei bulloni
for b, (x_pos, y_pos) in Bulloni.items():
    x1_cyl, y1_cyl, z1_cyl = create_cylinder(x_pos, y_pos, -tp-30, 0, db / 2)
    fig.add_trace(
        go.Mesh3d(
            x=x1_cyl, y=y1_cyl, z=z1_cyl,
            color='yellow', opacity=0.99, name=f"Bullone {b}"
        )
    )

    x2_cyl, y2_cyl, z2_cyl = create_cylinder(x_pos, y_pos, tp, 0, db / 2)
    fig.add_trace(
        go.Mesh3d(
            x=x2_cyl, y=y2_cyl, z=z2_cyl,
            color='yellow', opacity=0.5, name=f"Bullone {b}"
        )
    )
    #st.write(x2_cyl)
    for i in range(0, len(x1_cyl)-1):
        fig.add_trace(
            go.Mesh3d(
                x=[x1_cyl[i], x2_cyl[i], x2_cyl[i+1], x1_cyl[i+1]], 
                y=[y1_cyl[i], y2_cyl[i], y2_cyl[i+1], y1_cyl[i+1]],
                z=[z1_cyl[i], z2_cyl[i], z2_cyl[i+1], z1_cyl[i+1]],
                color='yellow', opacity=0.5, name=f"Bullone {b}",
        i = [0,1],
        j = [1,2],
        k = [3,3],
        ))

    ##DADO 1
    x3_cyl, y3_cyl, z3_cyl = create_cylinder(x_pos, y_pos, tp, 0, dd/2, 7)
    fig.add_trace(
        go.Mesh3d(
            x=x3_cyl, y=y3_cyl, z=z3_cyl,
            color='gray', opacity=0.5, name=f"Bullone {b}"
        )
    )

    x4_cyl, y4_cyl, z4_cyl = create_cylinder(x_pos, y_pos, tp+hd, 0, dd/2, 7)
    fig.add_trace(
        go.Mesh3d(
            x=x4_cyl, y=y4_cyl, z=z4_cyl,
            color='gray', opacity=0.5, name=f"Bullone {b}"
        )
    )
    #st.write(x2_cyl)
    for i in range(0, len(x3_cyl)-1):
        fig.add_trace(
            go.Mesh3d(
                x=[x3_cyl[i], x4_cyl[i], x4_cyl[i+1], x3_cyl[i+1]], 
                y=[y3_cyl[i], y4_cyl[i], y4_cyl[i+1], y3_cyl[i+1]],
                z=[z3_cyl[i], z4_cyl[i], z4_cyl[i+1], z3_cyl[i+1]],
                color='yellow', opacity=0.5, name=f"Bullone {b}",
        i = [0,1],
        j = [1,2],
        k = [3,3],
        ))

    ##DADO 2
    x3_cyl, y3_cyl, z3_cyl = create_cylinder(x_pos, y_pos, 0, 0, dd/2, 7)
    fig.add_trace(
        go.Mesh3d(
            x=x3_cyl, y=y3_cyl, z=z3_cyl,
            color='gray', opacity=0.5, name=f"Bullone {b}"
        )
    )

    x4_cyl, y4_cyl, z4_cyl = create_cylinder(x_pos, y_pos, -hd, 0, dd/2, 7)
    fig.add_trace(
        go.Mesh3d(
            x=x4_cyl, y=y4_cyl, z=z4_cyl,
            color='gray', opacity=0.5, name=f"Bullone {b}"
        )
    )
    #st.write(x2_cyl)
    for i in range(0, len(x3_cyl)-1):
        fig.add_trace(
            go.Mesh3d(
                x=[x3_cyl[i], x4_cyl[i], x4_cyl[i+1], x3_cyl[i+1]], 
                y=[y3_cyl[i], y4_cyl[i], y4_cyl[i+1], y3_cyl[i+1]],
                z=[z3_cyl[i], z4_cyl[i], z4_cyl[i+1], z3_cyl[i+1]],
                color='gray', opacity=0.5, name=f"Bullone {b}",
        i = [0,1],
        j = [1,2],
        k = [3,3],
        ))
    



# Aggiunta della freccia di forza
fig.add_trace(
    go.Scatter3d(
        x=[lpx + 10, lpx + 50],
        y=[lpy / 2, lpy / 2],
        z=[tp / 2, tp / 2],
        mode='lines',
        line=dict(color='rgb(255,51,0)', width=5),
        name="Freccia Forza"
    )
)
fig.add_trace(
    go.Cone(
        x=[lpx + 50], y=[lpy / 2], z=[tp / 2],
        u=[1], v=[0], w=[0],
        sizemode="absolute", sizeref=10,
        anchor="tip", colorscale=[[0, 'rgb(255,51,0)'], [1, 'rgb(255,51,0)']],
        name="Testa Freccia"
    )
)

# Configurazione del layout
fig.update_layout(
    scene=dict(
        xaxis=dict(title="Lunghezza (mm)", range=[-10, lpx + 60]),
        yaxis=dict(title="Larghezza (mm)", range=[-10, lpy + 10]),
        zaxis=dict(title="Spessore (mm)", range=[-10, tp + 10]),
        aspectmode='data',
    ),
    title="Piastra con Bulloni e Forza Applicata",
    width=800,
    height=600,
)

        
#---------------------------------------------------------------------------#        
#CALCOLO DELLE SOLLECITAZIONI

Lj = lpx - 2*e1

beta_lf = 1 - (Lj-15*db)/(200*db)

Ved_slu = V_slu/(nb*npt)
Ved_sle = V_sle/(nb*npt)

Fted_slu = T_slu/nb
Fted_sle = T_slu/nb

#forsa sulla piastra
Np = nf*Ved_slu

#---------------------------------------------------------------------------#
#CALCOLO DELLE RESISTENZE
Ares = selected_propB_value["Ares"]
fub = selected_classeB_value["fub"]

if tp > 40 and tp <= 80:
    fyk = selected_acciaio_value["fyk_40"]
    ftk = selected_acciaio_value["ftk_40"]
elif tp <= 40:
    fyk = selected_acciaio_value["fyk"]
    ftk = selected_acciaio_value["ftk"]

## Calcolo resistenza a taglio del singolo bullone
if selected_classeB_key == "4.6" or selected_classeB_key == "5.6" or selected_classeB_key == "8.8":
    Fvrd = (0.6*fub*Ares)/(gm2*1000)
if selected_classeB_key == "6.8" or selected_classeB_key == "10.9":
    Fvrd = (0.5*fub*Ares)/(gm2*1000)
    
##Calcolo resistenza a rifollamento
alpha_est = min(e1/(3*d0), fub/ftk, 1) #per bulloni di bordo in direzione del carico applicato
alpha_int = min(p1/(3*d0)-0.25, fub/ftk, 1) #per bulloni interni in direzione del carico applicato   
k_est = min(2.8*e2/(d0)-1.7, 2.5) #per bulloni di bordo in direzione perpendicolare al carico applicato
k_int = min(1.4*p2/(d0)-1.7, 2.5) #per bulloni interni in direzione perpendicolare al carico applicato   
Fbrd_est = k_est*alpha_est*ftk*db*tp/(gm2*1000) # resistenza per bullone esterno
Fbrd_int = k_int*alpha_int*ftk*db*tp/(gm2*1000) # resistenza per bullone interno

##Calcolo resistenza instabilità
res = 9/np.sqrt(235/fyk)

##Calcolo resistenza a trazione
Ftrd = 0.9*fub*Ares/(gm2*1000)

##Calcolo resistenza a punzonamento
dm = d0*1.0
Bprd = 0.6*np.pi*dm*tp*ftk/(gm2*1000)

##Calcolo resistenza ad attrito
selected_mu = st.sidebar.selectbox("Seleziona coeff. mu:", options=list(mu.keys()))
mu = mu[selected_mu]
Fpc = 0.7*Ares*fub/1000
Fpcd = Fpc/gm7_2
if T_slu == 0:
    Fsrd_slu = npt*mu*Fpcd/gm3_slu
else:
    Fsrd_slu = npt*mu*(Fpcd-0.8*Fted_slu)/gm3_slu

if T_sle == 0:
    Fsrd_sle = npt*mu*Fpcd/gm3_sle
else:
    Fsrd_sle = npt*mu*(Fpcd-0.8*Fted_sle)/gm3_slu
    


#resistenza sezione lorda della piastra
Rpl = lpy*tp*fyk/(gm0*1000) 
#resistenza sezione netta della piastra
Rpn = (lpy-d0*nf)*tp*fyk/(gm0*1000) 

#resistenza block shear
Lv = Lj+e1
if nf == 1:
    Lt = e1
else:
    Lt = p2

Anv = tp*(Lv-(d0*(nx-0.5)))
Ant = tp*(Lt- d0*nf*0.5 )
Agv = Lv*tp    
Veff1Rd = (ftk*Ant/gm2 + fyk*Anv/(np.sqrt(3)*gm0))/1000


    
    




# Configura layout del grafico
fig.update_layout(
    width=800,
    height=400,
    xaxis=dict(range=[-10, lpx + 10], title="Lunghezza (mm)", scaleanchor="y"),
    yaxis=dict(range=[-10, lpy + 10], title="Larghezza (mm)"),
    showlegend=False,
    title="Vista della Giunzione",
)


st.plotly_chart(fig)


st.subheader("Sintesi delle verifiche")
# Creiamo un dizionario con i dati delle verifiche
data = {
    "Verifica": [
        "Taglio del Bullone",
        "Rifollamento (Bullone di bordo)",
        "Rifollamento (Bullone interno)",
        "Instabilità del Piatto",
        "Trazione",
        "Punzonamento",
        "Taglio-Trazione",
        "Attrito (SLU)",
        "Attrito (SLE)",
        "sez. netta piastra",
        "block shear"
    ],
    
    # "Formula": [
    #     r"$F_{v,Rd} = \frac{\alpha_{s} \cdot A_{res} \cdot f_{tbk}}{\gamma_{M2}}$",
    #     r"$F_{b,Rd} = \frac{k \cdot \alpha \cdot f_{tk} \cdot d \cdot t}{\gamma_{M2}}$",
    #     r"$F_{b,Rd} = \frac{k \cdot \alpha \cdot f_{tk} \cdot d \cdot t}{\gamma_{M2}}$",
    #     r"$\frac{p_1}{t} ≤ 9\cdot \left(\frac{235}{f_{yk}}\right)^2$",
    #     r"$F_{t,Rd} = \frac{0.9 \cdot f_{tbk} \cdot A_{res}}{\gamma_{M2}}$",
    #     r"$B_{p,Rd} = \frac{0.6 \cdot \pi \cdot d_{m} \cdot t_{p} \cdot f_{tk}}{\gamma_{M2}}$",
    #     r"$\frac{F_{v,Ed}}{F_{v,Rd}} + \frac{F_{t,Ed}}{1.4F_{t,Rd}}$",
    #     r"$F_{s,Rd} = \frac{n \cdot \mu (F_{p,Cd} - 0.8F_{t,Ed})}{\gamma_{M2}}$",
    #     r"$F_{s,Rd} = \frac{n \cdot \mu (F_{p,Cd} - 0.8F_{t,Ed})}{\gamma_{M2}}$"
    # ],

    "D/C": [
        Ved_slu/Fvrd,
        Ved_slu/Fbrd_est,
        Ved_slu/Fbrd_int,
        (p1/tp)/(9/((235/fyk)**(0.5))),
        Fted_slu/Ftrd,
        Fted_slu/Bprd,
        (Ved_slu/Fvrd + Fted_slu/(1.4*Ftrd)),
        Ved_slu/Fsrd_slu,
        Ved_sle/Fsrd_sle,
        Np/Rpn,
        V_slu/Veff1Rd,
    ],
    "Esito": [
        "✅" if Ved_slu<= Fvrd else "❌",
        "✅" if Ved_slu<= Fbrd_est else "❌",
        "✅" if Ved_slu<= Fbrd_int else "❌",
        "✅" if p1/tp <= res else "❌ Instabilità",
        "✅" if Fted_slu/Ftrd <= 1.0 else "❌",
        "✅" if Fted_slu/Bprd <= 1.0 else "❌",
        "✅" if Ved_slu/Fvrd + Fted_slu/(1.4*Ftrd) <= 1 else "❌",
        "✅" if Ved_slu<= Fsrd_slu else "❌",
        "✅" if Ved_sle<= Fsrd_sle else "❌",
        "✅" if Np <= Rpn else "❌", #area netta
        "✅" if V_slu <= Veff1Rd else "❌", #block shear
    ]
}

# Creiamo un DataFrame con i dati
df = pd.DataFrame(data)

# Mostriamo la tabella
st.table(df)

# Output dei risultati
st.subheader("Risultati del Calcolo")
st.write(f"lunghezza della piastra: {lpx:.2f} mm")
st.write(f"larghezza della piastra: {lpy:.2f} mm")
st.write(f"forza sul bullone: {Ved_slu:.2f} KN")

st.subheader("Verifiche distanze ed interassi")
st.markdown(f" e1 = {e1:.2f} mm {'≤' if e1<= e1_min else '>'} e1_min = {e1_min:.2f} mm    {'✅ Verifica' if e1>= e1_min else '❌ Non Verifica'}")
st.markdown(f" p1 = {p1:.2f} mm {'≤' if p1<= p1_min else '>'} p1_min = {p1_min:.2f} mm    {'✅ Verifica' if p1>= p1_min else '❌ Non Verifica'}")
st.markdown(f" e2 = {e2:.2f} mm {'≤' if e2<= e2_min else '>'} e2_min = {e2_min:.2f} mm    {'✅ Verifica' if e2>= e2_min else '❌ Non Verifica'}")
st.markdown(f" p2 = {p2:.2f} mm {'≤' if p2<= p2_min else '>'} p2_min = {p2_min:.2f} mm    {'✅ Verifica' if p2>= p2_min else '❌ Non Verifica'}")

if Lj >= 15*d0:
    Fvrd = Fvrd*beta_lf


st.subheader("Verifica a taglio del bullone")

st.markdown("La resistenza al taglio di un bullone si calcola come:")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"F_{v,Rd} = \frac{\alpha_{s} \cdot A_{res} \cdot f_{tbk}}{\gamma_{M2}}")
# st.markdown("""
# - \ ftbk: Resistenza ultima del materiale del bullone.  
# - \ Ares: area resistente della parte filettata del bullone. 
# - \ alpha: 0.6 se il bullone è di classe 4.6, 5.6 e 8.8 oppure 0.5 se di classe 6.8 e 10.9.  
# - \( \gamma_{M2} \): Fattore parziale di sicurezza per i bulloni.
# """)
st.markdown(f" Ved = {Ved_slu:.2f} KN {'≤' if Ved_slu<= Fvrd else '>'} Fvrd = {Fvrd:.2f} KN ")
st.markdown(f" D/C = {Ved_slu/Fvrd:.2f} {'✅ Verifica' if Ved_slu<= Fvrd else '❌ Non Verifica'}")


st.subheader("Verifica a rifollamento")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"F_{b,Rd} = \frac{k \cdot \alpha \cdot f_{tk} \cdot d \cdot t}{\gamma_{M2}}")

st.markdown("**Verifica bullone di bordo**")
st.markdown(f" Ved = {Ved_slu:.2f} KN {'≤' if Ved_slu<= Fbrd_est else '>'} Fvrd = {Fbrd_est:.2f} KN ")
st.markdown(f" D/C = {Ved_slu/Fbrd_est:.2f} {'✅ Verifica' if Ved_slu<= Fbrd_est else '❌ Non Verifica'}")


st.markdown("**Verifica bullone interno**")
st.markdown(f" Ved = {Ved_slu:.2f} KN {'≤' if Ved_slu<= Fbrd_int else '>'} Fvrd = {Fbrd_int:.2f} KN ")
st.markdown(f" D/C = {Ved_slu/Fbrd_int:.2f} {'✅ Verifica' if Ved_slu<= Fbrd_int else '❌ Non Verifica'}")

st.subheader("Verifica instabilità del piatto")
st.markdown("La verifica è stata condotta rispettando la seguente equazione")
st.latex(r"\frac{p_1}{t} ≤ 9\cdot (235/f_{yk})^2")
st.markdown(f" p1/t = {p1/tp:.2f} {'<' if p1/tp < res else '>'} [9/(235/fyk)^0.5] = {res:.2f} {'✅ Verifica' if p1/tp<= res else '❌ deve essere considerata instabilità (l=0.6p1)'}")

st.subheader("Verifica trazione")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"F_{t,Rd} = \frac{0.9 \cdot f_{tbk} \cdot A_{res}}{\gamma_{M2}}")

st.markdown(f" Ted = {Fted_slu:.2f} KN {'≤' if Fted_slu<= Ftrd else '>'} Fsrd = {Ftrd:.2f} KN ")
st.markdown(f" D/C = {Fted_slu/Ftrd:.2f} {'✅ Verifica' if Fted_slu/Ftrd<= 1.0 else '❌ Non Verifica'}")

st.subheader("Verifica punzonamento")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"B_{p,Rd} = \frac{0.6 \cdot \pi \cdot d_{m} \cdot t_{p} \cdot f_{tk}}{\gamma_{M2}}")
st.markdown(f" Ted = {Fted_slu:.2f} KN {'≤' if Fted_slu <= Bprd else '>'} Bprd = {Bprd:.2f} KN ")
st.markdown(f" D/C = {Fted_slu/Bprd:.2f} {'✅ Verifica' if Fted_slu/Bprd<= 1.0 else '❌ Non Verifica'}")

st.subheader("Verifica taglio trazione")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"\frac{F_{v,Ed}}{F_{v,Rd}} + \frac{F_{t,Ed}}{1.4F_{t,Rd}}")
st.markdown(f" Ved/Vrd + Fted/(1.4Ftrd) = {Ved_slu/Fvrd + Fted_slu/(1.4*Ftrd):.2f} {'≤' if Ved_slu/Fvrd + Fted_slu/(1.4*Ftrd)<= 1 else '>'} 1 {'✅ Verifica' if Ved_slu/Fvrd + Fted_slu/(1.4*Ftrd) <= 1 else '❌ Non Verifica'}")

st.subheader("Verifica ad attrito")
st.markdown("La resistenza è valutata secondo $NTC18 4.2.8.1.1:")
st.latex(r"F_{s,Rd} = \frac{n \cdot \mu (F_{p,Cd} - 0.8F_{t,Ed})}{\gamma_{M2}}")
st.markdown("**Verifica ad attrito allo SLU**")
st.markdown(f" Ved = {Ved_slu:.2f} KN {'≤' if Ved_slu<= Fsrd_slu else '>'} Fsrd = {Fsrd_slu:.2f} KN ")
st.markdown(f" D/C = {Ved_slu/Fsrd_slu:.2f} {'✅ Verifica' if Ved_slu<= Fsrd_slu else '❌ Non Verifica'}")
st.markdown("**Verifica ad attrito allo SLE**")
st.markdown(f" Ved = {Ved_sle:.2f} KN {'≤' if Ved_sle<= Fsrd_sle else '>'} Fsrd = {Fsrd_sle:.2f} KN ")
st.markdown(f" D/C = {Ved_sle/Fsrd_sle:.2f} {'✅ Verifica' if Ved_sle<= Fsrd_sle else '❌ Non Verifica'}")


st.subheader("Verifica sezione netta della piastra")
st.markdown(f" Ned = {Np:.2f} KN {'≤' if Np <= Rpn else '>'} Nrd = {Rpn:.2f} KN ")
st.markdown(f" D/C = {Np/Rpn:.2f} {'✅ Verifica' if Np/Rpn<= 1.0 else '❌ Non Verifica'}")


st.subheader("Verifica block shear")
st.markdown(f" Ned = {V_slu:.2f} KN {'≤' if V_slu <= Veff1Rd else '>'} Nrd = {Veff1Rd:.2f} KN ")
st.markdown(f" D/C = {V_slu/Veff1Rd:.2f} {'✅ Verifica' if V_slu/Veff1Rd<= 1.0 else '❌ Non Verifica'}")

