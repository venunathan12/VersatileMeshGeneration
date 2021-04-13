#include <iostream>
#include <algorithm>
using namespace std;
#include <omp.h>

double MeshSize(double x, double y)
{
    double min = 1.0 / 32;
    return min + (x + y)*(x + y) / 8;
}

class Cell
{
public:
    double CenterX, CenterY, Size;
    Cell *TL, *TR, *BL, *BR;
    bool isParent, isActive;
    
    Cell(double x, double y, double s)
    {
        CenterX = x; CenterY = y; Size = s;
        TL = TR = BL = BR = nullptr;
        isParent = false; isActive = true;
    }
    Cell()
    {
        CenterX = 0.5; CenterY = 0.5; Size = 1.0;
        TL = TR = BL = BR = nullptr;
        isParent = false; isActive = true;
    }
};

#define MinDepth 2
#define MaxDepth 6

void SplitCell(Cell &ToSplit, Cell *Memory, int &Idx)
{
    ToSplit.isActive = false;
    Memory[Idx ++] = Cell(ToSplit.CenterX - ToSplit.Size / 4, ToSplit.CenterY - ToSplit.Size / 4, ToSplit.Size / 2);
    Memory[Idx ++] = Cell(ToSplit.CenterX + ToSplit.Size / 4, ToSplit.CenterY - ToSplit.Size / 4, ToSplit.Size / 2);
    Memory[Idx ++] = Cell(ToSplit.CenterX + ToSplit.Size / 4, ToSplit.CenterY + ToSplit.Size / 4, ToSplit.Size / 2);
    Memory[Idx ++] = Cell(ToSplit.CenterX - ToSplit.Size / 4, ToSplit.CenterY + ToSplit.Size / 4, ToSplit.Size / 2);
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

int main()
{
    Cell *Cells[MaxDepth + 1];
    int numCells[MaxDepth + 1];
    
    Cells[0] = new Cell[1];
    numCells[0] = 1;

    int Depth;
    for (Depth = 1; Depth <= MaxDepth; Depth ++)
    {
        Cells[Depth] = new Cell[4 * numCells[Depth - 1]];
        numCells[Depth] = 0;

        for (int i = 0; i < numCells[Depth - 1]; i++)
        {
            Cell &NowProcession = Cells[Depth - 1][i];
            
            if (Depth <= MinDepth)
            {
                SplitCell(NowProcession, Cells[Depth], numCells[Depth]);
            }
            else if (ShouldSplit(NowProcession))
            {
                SplitCell(NowProcession, Cells[Depth], numCells[Depth]);
            }
        }
    }

    for (Depth = 0; Depth <= MaxDepth; Depth ++)
        for (int i = 0; i < numCells[Depth]; i++)
            if (Cells[Depth][i].isActive)
                printf("%0.8lf %0.8lf %0.8lf\n", Cells[Depth][i].CenterX, Cells[Depth][i].CenterY, Cells[Depth][i].Size);

    return 0;
}