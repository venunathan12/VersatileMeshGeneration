#include <iostream>
#include <algorithm>
using namespace std;

#include <stdio.h>
#include <math.h>
#include <vector>
#include <stack>
#include <sys/time.h>

#include <omp.h>

#define MAX(l,r) (l>r?l:r)
#define MIN(l,r) (l<r?l:r)

#define NumThreads 4

#define MinDepth 2
#define MaxDepth 8

double MeshSize(double x, double y)
{
    double Min = 0, D = sqrt(x*x + y*y) - 0.5;
    return Min + MAX(D, -D) * 0.25;
}

class Point
{
public:
    double X, Y;

    Point(double x, double y)
    {
        X = x; Y = y;
    }
    Point()
    {
        X = Y = 0;
    }
};
bool operator < (Point L, Point R)
{
    if (L.Y != R.Y)
        return L.Y < R.Y;
    else
        return L.X < R.X;
}

class Cell
{
public:
    double CenterX, CenterY, Size;
    Cell *TL, *TR, *BL, *BR;
    Cell *L, *T, *R, *B;
    int IdC, IdBL, IdTL, IdTR, IdBR;
    bool isParent, isActive;
    
    Cell(double x, double y, double s)
    {
        CenterX = x; CenterY = y; Size = s;
        TL = TR = BL = BR = nullptr;
        isParent = false; isActive = true;
        L = T = R = B = nullptr;
        IdC = IdBL = IdTL = IdTR = IdBR = -1;
    }
    Cell()
    {
        CenterX = 0.5; CenterY = 0.5; Size = 1.0;
        TL = TR = BL = BR = nullptr;
        isParent = false; isActive = true;
        L = T = R = B = nullptr;
        IdC = IdBL = IdTL = IdTR = IdBR = -1;
    }
};

void SplitCell(Cell &ToSplit, Cell *Memory, int &Idx)
{
    ToSplit.isActive = false; ToSplit.isParent = true;
    ToSplit.BL = Memory + Idx; Memory[Idx ++] = Cell(ToSplit.CenterX - ToSplit.Size / 4, ToSplit.CenterY - ToSplit.Size / 4, ToSplit.Size / 2);
    ToSplit.BR = Memory + Idx; Memory[Idx ++] = Cell(ToSplit.CenterX + ToSplit.Size / 4, ToSplit.CenterY - ToSplit.Size / 4, ToSplit.Size / 2);
    ToSplit.TR = Memory + Idx; Memory[Idx ++] = Cell(ToSplit.CenterX + ToSplit.Size / 4, ToSplit.CenterY + ToSplit.Size / 4, ToSplit.Size / 2);
    ToSplit.TL = Memory + Idx; Memory[Idx ++] = Cell(ToSplit.CenterX - ToSplit.Size / 4, ToSplit.CenterY + ToSplit.Size / 4, ToSplit.Size / 2);

    ToSplit.BL -> R = ToSplit.BR; ToSplit.BR -> L = ToSplit.BL;
    ToSplit.TL -> R = ToSplit.TR; ToSplit.TR -> L = ToSplit.TL;
    ToSplit.BL -> T = ToSplit.TL; ToSplit.TL -> B = ToSplit.BL;
    ToSplit.BR -> T = ToSplit.TR; ToSplit.TR -> B = ToSplit.BR;
}
void ArrangeAttachments(Cell &ToSplit)
{
    if (ToSplit.L != nullptr)
        if (ToSplit.L -> isParent)
        {
            Cell &Neighbour = *(ToSplit.L);
            ToSplit.BL -> L = Neighbour.BR;
            ToSplit.TL -> L = Neighbour.TR;
        }
        else
        {
            ToSplit.BL -> L = ToSplit.TL -> L = ToSplit.L;
        }
    
    if (ToSplit.T != nullptr)
        if (ToSplit.T -> isParent)
        {
            Cell &Neighbour = *(ToSplit.T);
            ToSplit.TL -> T = Neighbour.BL;
            ToSplit.TR -> T = Neighbour.BR;
        }
        else
        {
            ToSplit.TL -> T = ToSplit.TR -> T = ToSplit.T;
        }
    
    if (ToSplit.R != nullptr)
        if (ToSplit.R -> isParent)
        {
            Cell &Neighbour = *(ToSplit.R);
            ToSplit.TR -> R = Neighbour.TL;
            ToSplit.BR -> R = Neighbour.BL;           
        }
        else
        {
            ToSplit.TR -> R = ToSplit.BR -> R = ToSplit.R;
        }
    
    if (ToSplit.B != nullptr)
        if (ToSplit.B -> isParent)
        {
            Cell &Neighbour = *(ToSplit.B);
            ToSplit.BR -> B = Neighbour.TR;
            ToSplit.BL -> B = Neighbour.TL;
        }
        else
        {
            ToSplit.BR -> B = ToSplit.BL -> B = ToSplit.B;
        }
}
bool ShouldSplit(Cell &ToDecide)
{
    double SizeEvaluations[5];
    SizeEvaluations[0] = MeshSize(ToDecide.CenterX, ToDecide.CenterY);
    SizeEvaluations[1] = MeshSize(ToDecide.CenterX - ToDecide.Size / 2, ToDecide.CenterY - ToDecide.Size / 2);
    SizeEvaluations[2] = MeshSize(ToDecide.CenterX + ToDecide.Size / 2, ToDecide.CenterY - ToDecide.Size / 2);
    SizeEvaluations[3] = MeshSize(ToDecide.CenterX + ToDecide.Size / 2, ToDecide.CenterY + ToDecide.Size / 2);
    SizeEvaluations[4] = MeshSize(ToDecide.CenterX - ToDecide.Size / 2, ToDecide.CenterY + ToDecide.Size / 2);

    if (* min_element(SizeEvaluations, SizeEvaluations + 5) < ToDecide.Size)
        return true;
    else
        return false;
}

class Triangle
{
public:
    int V1, V2, V3;

    Triangle()
    {
        V1 = V1 = V3 = -1;
    }

    Triangle(int v1, int v2, int v3)
    {
        V1 = v1, V2= v2, V3 = v3;
    }
};

vector <Point> Vertices;
vector <Triangle> Triangles;

int AddEntry(Point Entry, vector <Point> &Collection, vector <int> &IDs, int Id, int Sz)
{
    Collection.push_back(Entry);
    int PtId = Sz * Collection.size() + Id;
    IDs.push_back(PtId);
    return PtId;
}
int AddEntry(Triangle Entry, vector <Triangle> &Collection)
{
    Collection.push_back(Entry);
    return Collection.size();
}

void MarkCell(Cell &C, vector <Point> &Collection, vector <int> &IDs, int Id, int Sz)
{
    C.IdC = AddEntry(Point(C.CenterX, C.CenterY), Collection, IDs, Id, Sz);
    C.IdBL = AddEntry(Point(C.CenterX-C.Size/2, C.CenterY-C.Size/2), Collection, IDs, Id, Sz);

    bool SetTL = true;
    if (C.T != nullptr)
    {
        Cell &Nb = *(C.T);

        double NbStartX = Nb.CenterX - Nb.Size / 2;
        if (NbStartX == C.CenterX - C.Size / 2)
            SetTL = false;
    }
    if (SetTL && C.L == nullptr)
        C.IdTL = AddEntry(Point(C.CenterX-C.Size/2, C.CenterY+C.Size/2), Collection, IDs, Id, Sz);

    bool SetBR = true;
    if (C.R != nullptr)
    {
        Cell &Nb = *(C.R);

        double NbStartY = Nb.CenterY - Nb.Size / 2;
        if (NbStartY == C.CenterY - C.Size / 2)
            SetBR = false;
    }
    if (SetBR && C.B == nullptr)
        C.IdBR = AddEntry(Point(C.CenterX+C.Size/2, C.CenterY-C.Size/2), Collection, IDs, Id, Sz);

    bool SetTR = true;
    if (C.T != nullptr && C.R != nullptr)
    {
        Cell &NbT = *(C.T), &NbR = *(C.R);

        double NbTEndX = NbT.CenterX + NbT.Size / 2;
        double NbREndY = NbR.CenterY + NbR.Size / 2;
        if ((NbTEndX == C.CenterX + C.Size / 2) && (NbREndY == C.CenterY + C.Size / 2))
            SetTR = false;
    }
    if (SetTR)
        C.IdTR = AddEntry(Point(C.CenterX+C.Size/2, C.CenterY+C.Size/2), Collection, IDs, Id, Sz);
}
void CollectCorners(Cell &C)
{
    if (C.IdTR == -1)
    {
        Cell *S = C.T -> R;
        while (S -> isParent) S = S -> BL;
        C.IdTR = S -> IdBL;
    }

    if (C.IdTL == -1)
    {
        bool SetTL = true;
        if (C.T != nullptr)
        {
            Cell &Nb = *(C.T);

            double NbStartX = Nb.CenterX - Nb.Size / 2;
            if (NbStartX == C.CenterX - C.Size / 2)
                SetTL = false;
        }

        if (SetTL)
        {
            Cell *S = C.L;
            while (S -> isParent) S = S -> TR;
            C.IdTL = S -> IdTR;
        }
        else
        {
            Cell *S = C.T;
            while (S -> isParent) S = S -> BL;
            C.IdTL = S -> IdBL;
        }
    }

    if (C.IdBR == -1)
    {
        bool SetBR = true;
        if (C.R != nullptr)
        {
            Cell &Nb = *(C.R);

            double NbStartY = Nb.CenterY - Nb.Size / 2;
            if (NbStartY == C.CenterY - C.Size / 2)
                SetBR = false;
        }

        if (SetBR)
        {
            Cell *S = C.B;
            while (S -> isParent) S = S -> TR;
            C.IdBR = S -> IdTR;
        }
        else
        {
            Cell *S = C.R;
            while (S -> isParent) S = S -> BL;
            C.IdBR = S -> IdBL;
        }
    }
}
void MarkTris(Cell &C, vector <Triangle> &Collection)
{
    if (C.T != nullptr && C.T -> Size == C.Size)
    {
        stack <Cell *> TSt;
        TSt.push(C.T);

        while (! TSt.empty())
        {
            Cell *Z = TSt.top();
            TSt.pop();

            if (Z -> isParent)
                TSt.push(Z -> BR), TSt.push(Z -> BL);
            else
                AddEntry(Triangle(C.IdC, Z -> IdBR, Z -> IdBL), Collection);
        }
    }
    else
        AddEntry(Triangle(C.IdC, C.IdTR, C.IdTL), Collection);
    
    if (C.R != nullptr && C.R -> Size == C.Size)
    {
        stack <Cell *> RSt;
        RSt.push(C.R);

        while (! RSt.empty())
        {
            Cell *Z = RSt.top();
            RSt.pop();

            if (Z -> isParent)
                RSt.push(Z -> BL), RSt.push(Z -> TL);
            else
                AddEntry(Triangle(C.IdC, Z -> IdBL, Z -> IdTL), Collection);
        }
    }
    else
        AddEntry(Triangle(C.IdC, C.IdBR, C.IdTR), Collection);
    
    if (C.B != nullptr && C.B -> Size == C.Size)
    {
        stack <Cell *> BSt;
        BSt.push(C.B);

        while (! BSt.empty())
        {
            Cell *Z = BSt.top();
            BSt.pop();

            if (Z -> isParent)
                BSt.push(Z -> TL), BSt.push(Z -> TR);
            else
                AddEntry(Triangle(C.IdC, Z -> IdTL, Z -> IdTR), Collection);
        }
    }
    else
        AddEntry(Triangle(C.IdC, C.IdBL, C.IdBR), Collection);
    
    if (C.L != nullptr && C.L -> Size ==  C.Size)
    {
        stack <Cell *> LSt;
        LSt.push(C.L);

        while (! LSt.empty())
        {
            Cell *Z = LSt.top();
            LSt.pop();

            if (Z -> isParent)
                LSt.push(Z -> TR), LSt.push(Z -> BR);
            else
                AddEntry(Triangle(C.IdC, Z -> IdTR, Z -> IdBR), Collection);
        }
    }
    else
        AddEntry(Triangle(C.IdC, C.IdTL, C.IdBL), Collection);
}

long long GetTime()
{
    struct timeval T;
    gettimeofday(&T, NULL);
    return T.tv_usec + T.tv_sec * (long long) 1000000;
}

int main()
{
    long long Ts = GetTime();

    omp_set_num_threads(NumThreads);
    Cell BaseCell, **CellPtrs[MaxDepth + 1];
    int numCells[MaxDepth + 1], numCellsTotal = 0;
    int *localNums = new int[NumThreads + 1];
    Cell **allActiveCells;

    Point *FinalPoints; int *FinalIDs; Triangle *FinalTris;
    int numPtsTotal, numTrisTotal;

    CellPtrs[0] = new Cell*[1]; CellPtrs[0][0] = &BaseCell;
    numCells[0] = 1;

    #pragma omp parallel
    for (int Depth = 1; Depth <= MaxDepth; Depth ++)
    {
        int Tn = numCells[Depth - 1];
        int Id = omp_get_thread_num(), Sz = omp_get_num_threads();

        int ourSt = Id * (Tn / Sz) + MIN(Id, Tn % Sz);
        int ourEn = ourSt + Tn / Sz + (int)(Id < Tn % Sz);

        #pragma omp barrier

        Cell *ourMem = new Cell[4*(ourEn-ourSt)];
        int ourNumNew = 0;

        for (int i = ourSt; i < ourEn; i++)
        {
            Cell &NowProcession = *CellPtrs[Depth - 1][i];
            if (Depth <= MinDepth || ShouldSplit(NowProcession))
                SplitCell(NowProcession, ourMem, ourNumNew);
        }
        
        #pragma omp barrier
        
        for (int i = ourSt; i < ourEn; i++)
        {
            Cell &NowProcession = *CellPtrs[Depth - 1][i];
            if (NowProcession.isParent)
                ArrangeAttachments(NowProcession);
        }

        localNums[Id + 1] = ourNumNew;

        #pragma omp barrier

        #pragma omp master
        {
            localNums[0] = 0;
            for (int i = 1; i <= Sz; i++)
                localNums[i] += localNums[i-1];
            
            CellPtrs[Depth] = new Cell*[localNums[Sz]];
            numCells[Depth] = localNums[Sz];
        }

        #pragma omp barrier

        int nwSt = localNums[Id], nwEn = localNums[Id + 1];
        for (int i = nwSt; i < nwEn; i++)
            CellPtrs[Depth][i] = ourMem + i - nwSt;
    }

    #pragma omp parallel
    {
        int Id = omp_get_thread_num(), Sz = omp_get_num_threads();

        vector <Cell *> ourActiveCells;

        for (int Depth = 0; Depth <= MaxDepth; Depth ++)
        {
            #pragma omp for
            for (int i = 0; i < numCells[Depth]; i++)
                if (CellPtrs[Depth][i] -> isActive)
                    ourActiveCells.push_back(CellPtrs[Depth][i]);
        }

        localNums[Id + 1] = ourActiveCells.size();

        #pragma omp barrier

        #pragma omp master
        {
            localNums[0] = 0;
            for (int i = 1; i <= Sz; i++)
                localNums[i] += localNums[i-1];
            
            numCellsTotal = localNums[Sz];
            allActiveCells = new Cell*[numCellsTotal];
        }

        #pragma omp barrier

        int nwSt = localNums[Id], nwEn = localNums[Id + 1];
        for (int i = nwSt; i < nwEn; i++)
            allActiveCells[i] = ourActiveCells[i - nwSt];
        
        #pragma omp barrier

        vector <Point> ourPoints; vector <int> ourIDs;
        vector <Triangle> ourTris;

        #pragma omp for
        for (int i = 0; i < numCellsTotal; i++)
            MarkCell(*allActiveCells[i], ourPoints, ourIDs, Id, Sz);
        
        #pragma omp for
        for (int i = 0; i < numCellsTotal; i++)
            CollectCorners(*allActiveCells[i]);

        #pragma omp for
        for (int i = 0; i < numCellsTotal; i++)
            MarkTris(*allActiveCells[i], ourTris);
        
        localNums[Id + 1] = ourPoints.size();

        #pragma omp barrier

        #pragma omp master
        {
            localNums[0] = 0;
            for (int i = 1; i <= Sz; i++)
                localNums[i] += localNums[i-1];
            
            numPtsTotal = localNums[Sz];
            FinalPoints = new Point[numPtsTotal];
            FinalIDs = new int[numPtsTotal];
        }

        #pragma omp barrier

        nwSt = localNums[Id], nwEn = localNums[Id + 1];
        for (int i = nwSt; i < nwEn; i++)
            FinalPoints[i] = ourPoints[i - nwSt], FinalIDs[i] = ourIDs[i - nwSt];
        
        #pragma omp barrier

        localNums[Id + 1] = ourTris.size();

        #pragma omp barrier

        #pragma omp master
        {
            localNums[0] = 0;
            for (int i = 1; i <= Sz; i++)
                localNums[i] += localNums[i-1];
            
            numTrisTotal = localNums[Sz];
            FinalTris = new Triangle[numTrisTotal];
        }

        #pragma omp barrier

        nwSt = localNums[Id], nwEn = localNums[Id + 1];
        for (int i = nwSt; i < nwEn; i++)
            FinalTris[i] = ourTris[i - nwSt];
    }

    long long Tc = GetTime() - Ts;

    cout << "$MeshFormat" << endl;
    cout << "2.2 0 8" << endl;
    cout << "$EndMeshFormat" << endl << endl;

    cout << "$Nodes" << endl;
    cout << numPtsTotal << endl;
    for (int i = 0; i < numPtsTotal; i++)
        printf("%d %0.8lf %0.8lf %lf\n", FinalIDs[i], FinalPoints[i].X, FinalPoints[i].Y, 0.0);
    cout << "$EndNodes" << endl << endl;

    cout << "$Elements" << endl;
    cout << numTrisTotal << endl;
    for (int i = 0; i < numTrisTotal; i++)
        printf("%d 2 2 0 1 %d %d %d\n", i+1, FinalTris[i].V1, FinalTris[i].V2,FinalTris[i].V3);
    cout << "$EndElements" << endl;

    long long Tf = GetTime() - Ts;

    cout << endl;
    cout << "Computation completed in : " << Tc << " microseconds" << endl;
    cout << "Execution completed in : " << Tf << " microseconds" << endl;

    return 0;
}