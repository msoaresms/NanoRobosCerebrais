#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include <cstring>
using namespace std;

const bool DOENTE = true;
const bool SAUDAVEL = false;
const int NULO = -1;
const int INFINITO = numeric_limits<int>::max();
//----------------------------------------------------------------------------------------------------------------------
class UnionFind {
private:
    vector<int> conjuntos;
    unsigned int tamanho;

public:
    UnionFind(unsigned int pTamanho) {
        this->tamanho = pTamanho;
        this->conjuntos.resize(tamanho+1);
        for (unsigned int i = 0; i < this->conjuntos.size(); i++){
            this->conjuntos[i] = i;
        }
    }

    int procurar(int);
    bool conectados(int, int);
    void uniao(int, int);

};
//----------------------------------------------------------------------------------------------------------------------
class Vertice {
private:
    float distancia;
    int predecessor;
    int numVertice;

public:
    Vertice(){}
    Vertice(float pDistacia, int pPredecessor, int pNumVertice) {
        this->distancia = pDistacia;
        this->predecessor = pPredecessor;
        this->numVertice = pNumVertice;
    }

    float getDistancia() { return this->distancia; }
    void setDistancia(float pDistancia) { this->distancia = pDistancia; }
    int getPredecessor() { return this->predecessor; };
    void setPredecessor(int pPredecessor) { this->predecessor = pPredecessor; }
    int getNumVertice() { return this->numVertice; }
    void setNumVertice(int pNumVertice) { this->numVertice = pNumVertice; }
};
//----------------------------------------------------------------------------------------------------------------------
class Heap {
private:
    int tamanhoHeap, tamArray;

    int esquerda(int i) { return i * 2; }
    int direita(int i) { return ((i * 2) + 1); }
    int pai(int i) { return floor(i / 2); }
    void troca(Vertice *, int, int);

public:
    int getTamanhoHeap();
    void heapfica(Vertice *, int);
    void constroiHeap(Vertice *, int);
    Vertice extrairMenor(Vertice *);
    void alterarChave(Vertice *, int, float);
};
//----------------------------------------------------------------------------------------------------------------------
class Aresta {
private:
    int u, v;
    float peso;

public:
    Aresta() {}
    Aresta(int pU, int pV, float pPeso) {
        this->u = pU;
        this->v = pV;
        this->peso = pPeso;
    }

    int getU() { return this->u; }
    int getV() { return this->v; }
    float getPeso() { return this->peso; }
};
//----------------------------------------------------------------------------------------------------------------------
template <class T>
class Grafo {
private:
    vector<T> itens;
    vector<float> *matriz;
    int ordem, tamanho;

public:
    Grafo() { Grafo(0); }
    Grafo(int pOrdem) {
        this-> ordem = pOrdem;
        this->tamanho = 0;
        inicializar();
    }

    void inicializar();
    void inserirAresta(int, int, float);
    void inserirItem(T, int);
    T getItem(int);

    vector<float> *getMatriz() { return this->matriz; }
    int getOrdem() { return this->ordem; }
    int getTamanho() { return this->tamanho; }
};
//----------------------------------------------------------------------------------------------------------------------
int UnionFind::procurar(int p) {
    while (p != this->conjuntos[p]) { p = this->conjuntos[p]; }
    return p;
}

bool UnionFind::conectados(int p, int q) {
    return (this->procurar(p) == this->procurar(q));
}

void UnionFind::uniao(int p, int q) {
    int raizP = this->procurar(p);
    int raizQ = this->procurar(q);

    if (raizP == raizQ) { return; }

    this->conjuntos[raizP] = raizQ;
    this->tamanho--;
}
//----------------------------------------------------------------------------------------------------------------------
int Heap::getTamanhoHeap() { return this->tamanhoHeap; }

void Heap::troca(Vertice *v, int i, int j) {
    Vertice vAux = v[i];
    v[i] = v[j];
    v[j] = vAux;
}

void Heap::heapfica(Vertice *v, int i) {
    int esq = esquerda(i);
    int dir = direita(i);
    int menor;

    if (esq <= this->tamanhoHeap && (v[esq].getDistancia() < v[i].getDistancia())) {
        menor = esq;
    } else {
        menor = i;
    }

    if (dir <= this->tamanhoHeap && (v[dir].getDistancia() < v[menor].getDistancia())) {
        menor = dir;
    }

    if (menor != i) {
        this->troca(v, i, menor);
        this->heapfica(v, menor);
    }
}

void Heap::constroiHeap(Vertice *v, int pTamArray) {
    this->tamArray = pTamArray;
    this->tamanhoHeap = this->tamArray;

    for (int i = floor(this->tamArray/2); i > 0; i--) {
        this->heapfica(v, i);
    }
}

Vertice Heap::extrairMenor(Vertice *v) {
    Vertice menor = v[1];
    v[1] = v[this->tamanhoHeap];
    this->tamanhoHeap--;
    this->heapfica(v, 1);
    return menor;
}

void Heap::alterarChave(Vertice *v, int vert, float peso) {
    v[vert].setDistancia(peso);

    while(vert > 1 && v[pai(vert)].getDistancia() > v[vert].getDistancia()) {
        this->troca(v, vert, pai(vert));
        vert = pai(vert);
    }
}
//----------------------------------------------------------------------------------------------------------------------
template <typename T>
void Grafo<T>::inicializar() {
    this->itens.resize(this->ordem+1);
    this->matriz = new vector<float>[this->ordem+1];

    for (int i = 0; i <= this->ordem; i++) { matriz[i].assign(this->ordem+1, -1.0); }
}

template <typename T>
void Grafo<T>::inserirAresta(int u, int v, float w) {
    this->matriz[u][v] = w;
    this->matriz[v][u] = w;
    this->tamanho++;
}

template <typename T>
void Grafo<T>::inserirItem(T pItem, int pVertice) { this->itens[pVertice] = pItem; }

template <typename T>
T Grafo<T>::getItem(int n) { return this->itens[n]; }
//----------------------------------------------------------------------------------------------------------------------
int particao(vector<Aresta> &v, int p, int r){
    Aresta x = v[p];
    Aresta tmp = v[r+1];
    v[r+1] = x;
    int i = p;
    int j = r+1;

    while(true){
        do{
            i = i+1;
        }while((v[i].getPeso() < x.getPeso()) && (i <= r+1));
        do{
            j = j-1;
        }while ((v[j].getPeso() > x.getPeso()) && (j >= p));

        if (i < j){
            Aresta aux = v[i];
            v[i] = v[j];
            v[j] = aux;
        } else {
            Aresta aux = v[p];
            v[p] = v[j];
            v[j] = aux;
            v[r+1] = tmp;
            return j;
        }
    }
}

void quickSort(vector<Aresta> &v, int p, int r) {
    if (p < r){
        int q = particao(v, p, r);
        quickSort(v, p, q-1);
        quickSort(v, q+1, r);
    }
}

void quickSort(vector<Aresta> &v) {
    quickSort(v, 0, v.size()-1);
}

int posNoHeap(Vertice *v, int tam, int n) {
    for (int i = 1; i <= tam; i++) {
        if (n == v[i].getNumVertice()) {
            return i;
        }
    }
    return -1;
}

Grafo<Grafo<bool>> lerCerebro() {
    string entrada;
    vector<string> entradaSplit;
    int ordemCerebro, tamanhoCerebro = NULO;

    for (int i = 0; i < 2; i++) {
        cin >> entrada;
        entradaSplit.push_back(entrada);
    }
    ordemCerebro = stoi(entradaSplit[0]);
    tamanhoCerebro = stoi(entradaSplit[1]);

    Grafo<Grafo<bool>> cerebro(ordemCerebro);

    for (int i = 1; i <= tamanhoCerebro; i++) {
        entradaSplit.clear();
        for (int j = 0; j < 3; j++) {
            cin >> entrada;
            entradaSplit.push_back(entrada);
        }

        int u = stoi(entradaSplit[0]);
        int v = stoi(entradaSplit[1]);
        float w = stof(entradaSplit[2]);

        cerebro.inserirAresta(u, v, w);
    }
    return cerebro;
}

void lerBlocos(Grafo<Grafo<bool>> &cerebro) {
    string entrada;
    vector<string> entradaSplit;
    int ordemCerebro = cerebro.getOrdem();

    for (int i = 1; i <= ordemCerebro; i++) {
        int ordemBloco, tamanhoBloco = NULO;
        vector<int> doentes;
        int nDoentes = 0;

        entradaSplit.clear();
        for (int j = 0; j < 2; j++) {
            cin >> entrada;
            entradaSplit.push_back(entrada);
        }
        ordemBloco = stoi(entradaSplit[0]);
        tamanhoBloco = stoi(entradaSplit[1]);

        Grafo<bool> bloco(ordemBloco);

        cin >> entrada;
        nDoentes = stoi(entrada);

        if (nDoentes != 0) {
            entradaSplit.clear();
            for (int j = 0; j < nDoentes; j++) {
                cin >> entrada;
                entradaSplit.push_back(entrada);
                doentes.push_back(stoi(entradaSplit[0]));
            }
        }

        for (int j = 1; j <= tamanhoBloco; j++) {
            entradaSplit.clear();
            for (int k = 0; k < 3; k++) {
                cin >> entrada;
                entradaSplit.push_back(entrada);
            }

            int u = stoi(entradaSplit[0]);
            int v = stoi(entradaSplit[1]);
            float peso = stof(entradaSplit[2]);

            bloco.inserirAresta(u, v, peso);
        }

        for (int j = 1; j <= ordemBloco; j++) {
            bloco.inserirItem(false, j);
        }

        if (nDoentes != 0) {
            for (int x : doentes) {
                bloco.inserirItem(true, x);
            }
        }

        cerebro.inserirItem(bloco, i);

    }
}

vector<unsigned int> adjDeU(vector<float> aux) {
    vector<unsigned int> adjU;
    for (unsigned int i = 1; i < aux.size(); i++) {
        if (aux[i] >= 0.0) {
            adjU.push_back(i);
        }
    }
    return adjU;
}

template <typename T>
void Dijkstra(Grafo<T> grafo, vector<Vertice> &S,  int inicio, int fim) {
    Vertice Q[grafo.getOrdem()+1];
    vector<unsigned int> adjU;
    vector<float> *matriz = grafo.getMatriz();
    Heap prioridade;
    for (int i = 1; i <=grafo.getOrdem(); i++) {
        Q[i] = Vertice(INFINITO, NULO, i);
    }

    Q[inicio].setDistancia(0.0);
    prioridade.constroiHeap(Q, grafo.getOrdem());

    while (prioridade.getTamanhoHeap() > 0) {
        Vertice u = prioridade.extrairMenor(Q);
        S.push_back(u);

        if (u.getNumVertice() == fim) {
            return;
        }

        adjU = adjDeU(matriz[u.getNumVertice()]);
        for (int x : adjU) {
            //int pos = prioridade.posNoHeap(x);
            int pos = posNoHeap(Q, prioridade.getTamanhoHeap(), x);
            if (pos != -1 && Q[pos].getDistancia() > u.getDistancia() + (matriz[u.getNumVertice()][Q[pos].getNumVertice()])) {

                float distancia = u.getDistancia() + (matriz[u.getNumVertice()][Q[pos].getNumVertice()]);
                Q[pos].setPredecessor(u.getNumVertice());

                prioridade.alterarChave(Q, pos, distancia);
            }
        }
    }
}

template <typename T>
vector<Aresta> Kruskal(Grafo<T> grafo) {
    vector<Aresta> A;
    vector<Aresta> saida;
    UnionFind conjuntos(grafo.getOrdem());

    vector<float> *matriz = grafo.getMatriz();

    for (int i = 1; i <= grafo.getOrdem(); i++) {
        for (int j = 1; j <= grafo.getOrdem(); j++) {
            if (j > i && matriz[i][j] >= 0.0) {
                A.push_back(Aresta(i, j, matriz[i][j]));
            }
        }
    }

    quickSort(A);

    while (A.size() > 0) {
        Aresta aux = A[0];
        A.erase(A.begin()+0);

        if (!conjuntos.conectados(aux.getU(), aux.getV())) {
            saida.push_back(aux);
            conjuntos.uniao(aux.getU(), aux.getV());
        }
    }

    return saida;

}
//----------------------------------------------------------------------------------------------------------------------
int main() {
    string entrada;
    vector<string> entradaSplit;
    Grafo<Grafo<bool>> cerebro;
    int inicio, fim = NULO;

    cerebro = lerCerebro();

    for (int i = 0; i < 2; i++) {
        cin >> entrada;
        entradaSplit.push_back(entrada);
    }
    inicio = stoi(entradaSplit[0]);
    fim = stoi(entradaSplit[1]);

    lerBlocos(cerebro);
    vector<Vertice> resultD;
    Dijkstra(cerebro, resultD, inicio, fim);

    //------------------------------------------------------------------------------------------------------------------

    vector<int> blocosV;
    int aux = fim;
    while (aux != NULO) {
        for (Vertice x : resultD) {
            if (aux == x.getNumVertice()) {
                blocosV.push_back(x.getNumVertice());
                aux = x.getPredecessor();
                break;
            }
        }
    }
    vector<Grafo<bool>> blocos;
    for (int x : blocosV) {
        Grafo<bool> bloco = cerebro.getItem(x);
        for (int i = 1; i <= bloco.getOrdem(); i++) {
            bool aux = bloco.getItem(i);
            if (aux) {
                blocos.push_back(bloco);
                break;
            }
        }
    }

    float saida = 0.0;

    for (Grafo<bool> x : blocos) {
        vector<Aresta> resultKruskal = Kruskal(x);

        for (Aresta y : resultKruskal) {
            saida += y.getPeso();
        }
    }
    cout << saida << endl;
    return 0;
}