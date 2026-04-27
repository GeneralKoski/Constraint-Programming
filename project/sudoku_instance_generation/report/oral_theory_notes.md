# Note di Teoria per l'Orale

Riassunto rapido degli argomenti teorici del corso che possono essere chiesti durante la discussione del progetto Sudoku Instance Generation. Ogni sezione contiene la "lavagna mentale" minima da poter spiegare a voce.

---

## 1. Teorema di Régin e filtering di `alldifferent`

### Setup

Il vincolo `alldifferent(x₁, …, xₙ)` impone che tutte le variabili abbiano valori distinti. Una formulazione naive è `forall(i ≠ j) (xᵢ ≠ xⱼ)`, ma propaga male: ottiene solo bounds consistency, non arc consistency.

### Costruzione del grafo bipartito

Si costruisce il **value graph** `G = (V, E)`:
- nodi a sinistra: variabili `x₁, …, xₙ`;
- nodi a destra: valori `v₁, …, vₘ` nei domini;
- archi: `(xᵢ, v) ∈ E ⟺ v ∈ D(xᵢ)`.

Una **soluzione** del vincolo è un *matching* che copre tutte le variabili (perfect matching dal lato delle variabili).

### Berge (1970) e Régin (1994)

**Teorema (Berge)**: un matching `M` è massimo ⟺ non esiste un *augmenting path* rispetto a `M`.

**Conseguenza utile**: un arco `(xᵢ, v)` appartiene a *qualche* matching massimo ⟺ appartiene a `M` o a un *alternating cycle* o *alternating path* che parte da un nodo libero.

### Algoritmo di filtering (Régin, 1994)

1. trovare un matching massimo `M` (Hopcroft-Karp, `O(E·√V)`);
2. se `|M| < n`: il vincolo è **inconsistente** → fail;
3. orientare gli archi: archi di `M` da right→left, gli altri da left→right;
4. trovare le **componenti fortemente connesse** (SCC) e i nodi raggiungibili da nodi liberi;
5. un arco `(xᵢ, v)` può essere rimosso dal dominio sse non sta in `M` e non è in una SCC con un arco di `M`, e non è raggiungibile da un nodo right libero via alternating path.

**Complessità**: `O(E + V)` per il filtering (dopo aver pre-calcolato il matching). Permette di mantenere arc consistency in modo efficiente.

### Relevance per Sudoku

Ogni riga, colonna, blocco è un `alldifferent` su 9 variabili e 9 valori. Régin garantisce che ogni propagazione elimina dai domini *tutti* i valori incoerenti, non solo i bound. Su una griglia con poche celle compilate, questo accelera molto la propagazione iniziale e il branching.

---

## 2. NP-completezza del Sudoku generalizzato

### Statement

Il problema decisionale "data una griglia parziale Sudoku `n²×n²` (con blocchi `n×n`), esiste un completamento valido?" è **NP-completo**.

### Dimostrazione (sketch)

**NP**: dato un completamento, si verifica in tempo polinomiale che ogni riga, colonna e blocco contengano valori distinti.

**NP-hardness**: si riduce **Latin Square Completion** (NP-completo, Colbourn 1984) a Sudoku. Una `n²×n²` Sudoku partial-filled può codificare una `n²×n²` Latin square, perché i vincoli di blocco si possono "soddisfare banalmente" se aggiungiamo righe ausiliarie. Più preciso: ogni Latin square `n×n` può essere completato sse esiste un Sudoku `n²×n²` di appropriata struttura.

Yato e Seta (2003) hanno mostrato anche che il problema "ha esattamente una soluzione?" (USAT su Sudoku) è in **DP**, e contare le soluzioni è **#P-completo**.

### Implicazione pratica

Per Sudoku 9×9 fissato (n=3), il problema è banale (numero di griglie limitato), ma l'analisi asintotica giustifica perché serve un solver CP serio invece di un'enumeration brute-force scalabile.

---

## 3. #P-completezza del solution counting

### Definizione

Una funzione `f` è **#P** se esiste una macchina di Turing non-deterministica polinomiale `M` tale che `f(x) = #(percorsi accettanti di M su x)`.

`#P-completa` per riduzione parsimoniosa.

### Per Sudoku

"Quante completamenti validi ha questa griglia parziale?" è #P-completo (Yato 2003). Implicazione: non si conosce un algoritmo polinomiale per contare; al peggio si enumera.

### Cosa significa per il progetto

Per il check di unicità non serve contare *tutte* le soluzioni: basta sapere se ce ne sono `≥2`. Questo è un caso facile (basta fermarsi al secondo). Quindi nella pratica si fa una "approssimazione one-sided": il counting si limita a `n=2` invece di esaustivo. Tecnicamente: `#≥2` è in NP (basta esibire due soluzioni), non #P-completo.

Solve-and-block fa la stessa cosa con un'altra struttura: solve, blocking, solve. Anche questo è in NP.

---

## 4. Arc Consistency vs Bounds Consistency

### Definizioni

Sia `c(x, y)` un vincolo binario.

**Arc-consistent (AC)**: per ogni `a ∈ D(x)`, esiste `b ∈ D(y)` con `c(a, b)`. E viceversa.

**Bounds-consistent (BC)**: per ogni `a ∈ {min(D(x)), max(D(x))}`, esiste `b ∈ D(y)` con `c(a, b)`. E simmetricamente.

AC ⟹ BC, ma non viceversa.

### Algoritmi

- **AC-1** (Mackworth 1977): rivisita ripetutamente tutti i vincoli finché niente cambia. `O(end³)` con `n` variabili, `e` vincoli, `d` dim. dominio. Inefficace.
- **AC-3** (Mackworth 1977): mantiene una coda di archi da rivedere. `O(ed³)`. Standard nei solver moderni.
- **AC-4** (Mohr e Henderson 1986): tiene supporti pre-calcolati. `O(ed²)` ma overhead alto in pratica.

### Per Sudoku

Su `alldifferent`, AC è realizzato in modo *globale* dall'algoritmo di Régin (vedi §1), non come AC-3 sui singoli `≠`. La differenza è netta:

- formulazione disgiuntiva con `xᵢ ≠ xⱼ` ottiene solo BC su singoli vincoli (quando `D(xᵢ) = {a}`, rimuove `a` da `D(xⱼ)`);
- `alldifferent` come global ottiene AC: rimuove ogni valore non sopportato da un matching massimo.

In pratica su Sudoku la differenza è enorme nei branch iniziali, perché molte celle hanno dominio pieno `{1..9}` e le sole disuguaglianze a coppie non bastano a propagare cancellazioni globali.

---

## 5. Solve-and-block vs Solution counting (sintesi)

| Aspetto | Solve-and-block | Solution counting |
|---|---|---|
| Numero di ricerche | 2 (con re-start) | 1 (continua) |
| Vincoli aggiuntivi | sì (blocking constraint) | no |
| Esprimibile in MiniZinc | con un secondo modello | con flag `-a -n 2` |
| Dimostrare unicità | UNSAT della seconda ricerca | esaurire l'albero |
| Overhead per puzzle unici | doppio start | albero completo |
| Vantaggio | propagazione dalla blocking | riuso del search state |

Nel progetto: counting ~1.8× più veloce su Python; entrambi corretti.

---

## 6. Domande probabili e risposte rapide

**"Perché alldifferent è meglio di tante disuguaglianze?"**
→ AC vs BC. Régin filtra valori non sopportati da nessun matching massimo, le disuguaglianze a coppie no. Vedi §1.

**"Perché due metodi di unicità invece di uno?"**
→ Per *confrontarli* come richiesto dalla specifica. Il counting riusa search state, solve-and-block riparte. Sono entrambi corretti, differiscono in efficienza. Vedi §5.

**"Come gestisci il caso UNKNOWN da timeout?"**
→ Mai trattato come unique. Si rolledback la rimozione. Si registra nei log per quantificare. Vedi report §3.3.

**"Il minimo a 17 indizi vale per le tue strategie?"**
→ No, le strategie greedy senza backtracking sulla rimozione si fermano a 22-25. Il minimo a 17 (McGuire 2012) richiede algoritmi più sofisticati e ore di calcolo. Vedi report §6.2.

**"Qual è il ruolo del backend Python?"**
→ Validazione end-to-end senza dipendere da MiniZinc, e benchmark veloce (l'overhead di startup MiniZinc è ~400ms). Non sostituisce il modello CP, lo riproduce procedurally per il controllo di unicità. Vedi report §5.1.

**"Perché generare le griglie complete invece di usare solo Kaggle?"**
→ Per autonomia: il dataset Kaggle è un'opzione, ma il modello CP `sudoku_generate_full_grid.mzn` con `indomain_random` produce griglie distinte per ogni seed. Su 50 seed otteniamo 50 griglie uniche.

**"Cosa cambierebbe con search annotations diverse?"**
→ Con dom_w_deg si pesa la storia dei fallimenti, utile su istanze fortemente vincolate. Su Sudoku 9×9 il guadagno è marginale; su 16×16 sarebbe più evidente.

**"Quanto sono utili i vincoli ridondanti `sum=45`?"**
→ Marginale su Sudoku 9×9 perché alldifferent è già forte. Servono per dimostrazione pedagogica: i vincoli ridondanti rinforzano la propagazione lineare e aiutano in problemi con meno propagator nativo.
