using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

class DoubleSearch
{
    //***********************************************************************
    string graphstr = "";//ДЛЯ ЗАГРУЗКИ ИЗ ФАЙЛА
    int Nodes = 0; //узлы
    int Edges = 0; //дуги
    int Begin_Node = 0; //начальная вершина
    int K_Path = 0; //количество искомых путей
    int iter = 0;
    static double[] sum;//Для сравнения
    int[] arr;//ДЛЯ ЗАГРУЗКИ ИЗ ФАЙЛА
    Matrix<double> MatrixINF;//МАТРИЦА БЕСКОНЕЧНОСТИ
    Matrix<double> matrixx0;// МАТРИЦА (0,?,?)
    Matrix<double>[,] D_Matrix;//МАТРИЦА D
    Matrix<double>[,] L_Matrix;
    Matrix<double>[,] U_Matrix;
    Tuple<int, int, double> gr = null;// КОРТЕЖ УЗЕЛ 1 УЗЕЛ 2 ВЕС РЕБРА
    Dictionary<int, Tuple<int, int, double>> graph;//узел 1 узел 2 вес ребра
    Matrix<double>[] D_straight;//матрица векторов прямого поиска
    Matrix<double>[] D_straightCopy;//копия для вычислений
    Matrix<double>[] D_back;//матрица векторов обратного поиска
    Matrix<double>[] D_backCopy;//копия для вычислений
    int count;//счетчик шагов
    Task task;
    //***********************************************************************

    /// <summary>
    /// Алгоритм работы
    /// </summary>
    void ALGORITHM()
    {
        graph = new Dictionary<int, Tuple<int, int, double>>();
        LoadGraph();
        //
        MatrixINF = CreateMatrix.Dense(1, K_Path, (x, y) => double.PositiveInfinity);
        D_Matrix = new Matrix<double>[Nodes, Nodes];
        task = new Task(() => Initialize_D());
        task.Start();
        task.Wait();
        task = null;
        //Initialize_D();
        //
        matrixx0 = CreateMatrix.Dense(1, K_Path, (x, y) => double.PositiveInfinity);
        matrixx0[0, 0] = 0;
        arr = new int[3];
        L_Matrix = new Matrix<double>[Nodes, Nodes];
        U_Matrix = new Matrix<double>[Nodes, Nodes];
        //
        task = new Task(() => initialize_L_U());
        task.Start();
        task.Wait();
        task = null;
        //initialize_L_U();
        //
        OutputInitialData();
        Console.WriteLine("Матрица D:");
        Output_matr(D_Matrix);
        Console.WriteLine("Нижнетреугольная матрица L");
        Output_matr(L_Matrix);
        Console.WriteLine("Верхнетреугольная матрица U:");
        Output_matr(U_Matrix);
        Console.WriteLine("Алгоритм поиска:");

        D_straight = new Matrix<double>[Nodes];
        D_back = new Matrix<double>[Nodes];
        D_backCopy = new Matrix<double>[Nodes];
        D_straightCopy = new Matrix<double>[Nodes];

        Parallel.For(0, Nodes, i =>
        {
            D_straight[i] = MatrixINF;
            D_back[i] = MatrixINF;
            D_backCopy[i] = MatrixINF;
            D_straightCopy[i] = MatrixINF;
        });

        D_straight[0] = matrixx0;
        D_back[0] = matrixx0;
        int strok = 0;
        count = 0;
        do
        {
            ReverseLookup(strok);
            PrintBack();
            if (D_back[1].ToRowArrays()[0].ToArray().SequenceEqual(D_straight[1].ToRowArrays()[0].ToArray())
& D_back[2].ToRowArrays()[0].ToArray().SequenceEqual(D_straight[2].ToRowArrays()[0].ToArray())
& D_back[3].ToRowArrays()[0].ToArray().SequenceEqual(D_straight[3].ToRowArrays()[0].ToArray()) & count != 1)
                break;
            DirectSearch(strok);
            PrintStraight();
        }
        while (D_back != D_straight);

        PrintAnsw();
    }
    /// <summary>
    /// Вывести количество шагов и минимальные пути
    /// </summary>
    private void PrintAnsw()
    {
        Console.WriteLine("Количество шагов - " + count);
        Console.WriteLine("Ответ:");
        List<double> sum = new List<double>();
        for (int i = 0; i < D_back.Length; i++)
        {
            sum.Add(0);
            for (int j = 0; j < D_back[i].ToRowArrays()[0].ToArray().Length; j++)
            {
                sum[i] += Math.Abs(D_back[i].ToRowArrays()[0].ToArray()[j]);
            }
        }
        iter = 0;
        double min = sum[0];
        int[] number = new int[K_Path];
        while (iter != K_Path)
        {
            min = sum[0];
            for (int i = 0; i < sum.Count; ++i)
            {
                if (sum[i] < min)
                {
                    min = sum[i];
                    number[iter] = i;
                }
            }
            sum.Remove(sum[number[iter]]);
            iter++;
        }
        for (int i = 0; i < K_Path; i++)
        {
            Console.Write("[");
            for (int j = 0; j < D_back[i].ToRowArrays()[0].ToArray().Length; j++)
            {
                Console.Write(D_back[number[i]].ToRowArrays()[0].ToArray()[j] + " ");
            }
            Console.Write("]" + "\t");
        }
    }

    /// <summary>
    /// Вывести шаг обратного хода
    /// </summary>
    private void PrintBack()
    {
        Console.Write("Обратный:");
        for (int i = 0; i < Nodes; i++)
        {
            for (int j = 0; j < D_back[i].ToRowArrays()[0].ToArray().Length; j++)
            {
                Console.Write(D_back[i].ToRowArrays()[0].ToArray()[j] + " ");

            }
            Console.Write("\t");
        }
        Console.WriteLine();
    }

    /// <summary>
    /// Вывести шаг прямого хода
    /// </summary>
    private void PrintStraight()
    {
        Console.Write("Прямой:");
        for (int i = 0; i < Nodes; i++)
        {
            for (int j = 0; j < D_straight[i].ToRowArrays()[0].ToArray().Length; j++)
            {
                Console.Write(D_straight[i].ToRowArrays()[0].ToArray()[j] + " ");
            }
            Console.Write("\t");
        }
        Console.WriteLine();
    }

    /// <summary>
    /// Вывод начальных данных
    /// </summary>
    private void OutputInitialData()
    {
        Console.WriteLine("Количество узлов:" + Nodes + "\n" + "Колисчество дуг:" + Edges + "\n" + "Начальный узел из которого будет осуществляться поиск:" + Begin_Node + "\n" + "Количество искомых путей:" + K_Path);
        Console.WriteLine("Граф(вершина вершина вес):");
        foreach (var VARIABLE in graph.Values)
        {
            Console.Write(VARIABLE.Item1 + " " + VARIABLE.Item2 + " " + VARIABLE.Item3);
            Console.WriteLine();
        }
    }

    /// <summary>
    /// Вывести граф 
    /// </summary>
    private void Output_matr(Matrix<double>[,] D)
    {
        for (int i = 0; i < Nodes; i++)
        {
            for (int j = 0; j < Nodes; j++)
            {
                for (int k = 0; k < D_Matrix[i, j].ToArray().Length; k++)
                {
                    Console.Write(D[i, j].ToArray()[0, k] + " ");
                }

                Console.Write("\t");
            }

            Console.WriteLine();
        }
    }

    /// <summary>
    /// Загружаем граф из файла
    /// </summary>
    private void LoadGraph()
    {
        using (var streamreader = new StreamReader(@"C:\Users\Admin\Desktop\DoubleSearch_dorefactorit`1\DoubleSearch\graph.txt"))
        {
            Nodes = Int32.Parse(streamreader.ReadLine());
            Edges = Int32.Parse(streamreader.ReadLine());
            Begin_Node = Int32.Parse(streamreader.ReadLine());
            K_Path = Int32.Parse(streamreader.ReadLine());
            int EdgesCopy = Edges;
            while (EdgesCopy != 0)
            {
                graphstr = streamreader.ReadLine();
                arr = graphstr.Split(' ').Select(n => int.Parse(n)).ToArray();
                gr = new Tuple<int, int, double>(arr[0], arr[1], arr[2]);
                graph.Add(iter, gr);
                iter++;
                EdgesCopy--;
            }
        }
    }
    /// <summary>
    /// Заполняем матрицу D
    /// </summary>
    /// <param name="tuple"></param>
    private void Initialize_D()
    {
        matrixx0 = CreateMatrix.Dense(1, K_Path, (x, y) => double.PositiveInfinity);
        matrixx0[0, 0] = 0;
        Parallel.For(0, Nodes, i =>
        {
            for (int j = 0; j < Nodes; j++)
            {
                if (i == j)
                {
                    D_Matrix[i, j] = matrixx0;
                }
                if (D_Matrix[i, j] == null)
                    D_Matrix[i, j] = MatrixINF;
            }
        });
        matrixx0 = null;
        var tuple = graph.GetEnumerator();
        while (tuple.MoveNext() != false)
        {
            matrixx0 = CreateMatrix.Dense(1, K_Path, (x, y) => double.PositiveInfinity);
            matrixx0[0, 0] = tuple.Current.Value.Item3;
            if (D_Matrix[tuple.Current.Value.Item1 - 1, tuple.Current.Value.Item2 - 1].ToRowArrays()[0].ToArray()[0] != double.PositiveInfinity)
            {
                int i = 0;
                while (D_Matrix[tuple.Current.Value.Item1 - 1, tuple.Current.Value.Item2 - 1].ToRowArrays()[0].ToArray()[i] != double.PositiveInfinity)
                { i++; }
                D_Matrix[tuple.Current.Value.Item1 - 1, tuple.Current.Value.Item2 - 1][0, i] = tuple.Current.Value.Item3;
            }
            else
                D_Matrix[tuple.Current.Value.Item1 - 1, tuple.Current.Value.Item2 - 1] = matrixx0;
        }
    }
    /// <summary>
    /// Инициализация матриц L и U
    /// </summary>
    private void initialize_L_U()
    {
        Parallel.Invoke(
            () =>
            {
                for (int i = 0; i < Nodes; i++)
                {
                    for (int j = 0; j < Nodes; j++)
                    {
                        if (i < j)
                            U_Matrix[i, j] = D_Matrix[i, j];
                        else
                            U_Matrix[i, j] = MatrixINF;
                    }
                }
            },
            () =>
            {
                for (int i = 0; i < Nodes; i++)
                {
                    for (int j = 0; j < Nodes; j++)
                    {
                        if (i > j)
                            L_Matrix[i, j] = D_Matrix[i, j];
                        else
                            L_Matrix[i, j] = MatrixINF;
                    }
                }
            });
    }
    /// <summary>
    /// прямой поиск
    /// </summary>
    /// <param name="strok"></param>
    private async void DirectSearch(int strok)
    {
        Parallel.For(0, Nodes, i =>
        {
            D_straight[i] = D_back[i];
        });
        count++;
        strok = 0;
        for (int i = 0; i < Nodes; i++)
        {
            D_straightCopy[i] = MatrixINF;
        }

        for (int i = 0; i < Nodes; i++)
        {
            for (int j = 0; j < Nodes; j++)
            {
                if (j < i)
                {
                    if (U_Matrix[j, i] != MatrixINF)
                    {
                        await AsyncCalculationStraight(strok, i, j);
                        strok++;
                    }
                    else
                    {
                        strok++;
                        continue;
                    }
                }
            }
            strok = 0;
            D_straight[i] = await Add(D_straightCopy[i], D_back[i]);
        }
    }
    /// <summary>
    /// Асинхронный вспомогательный метод для прямого поиска
    /// </summary>
    /// <param name="strok"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <returns></returns>
    async private Task AsyncCalculationStraight(int strok, int i, int j)
    {
        D_straightCopy[i] = await Add(D_straightCopy[i], await Multipl(U_Matrix[strok, i], D_straight[j], K_Path));
    }

    /// <summary>
    /// Обратный поиск
    /// </summary>
    /// <param name="strok"></param>
    private async void ReverseLookup(int strok)
    {
        Parallel.For(0, Nodes, i =>
        {
            D_back[i] = D_straight[i];
        });
        strok = Nodes - 1;
        count++;
        for (int i = 0; i < Nodes; i++)
        {
            D_backCopy[i] = MatrixINF;
        }

        for (int i = Nodes - 1; i >= 0; i--)
        {
            for (int j = Nodes - 1; j >= 0; j--)
            {
                if (j > i)
                {
                    if (L_Matrix[j, i] != MatrixINF)
                    {
                        await AsyncCalcualtion1Back(strok, i, j);
                        strok--;
                    }
                    else
                    {
                        strok--;
                        continue;
                    }
                }
            }
            strok = Nodes - 1;
            D_back[i] = await Add(D_back[i], D_backCopy[i]);
        }
    }
    /// <summary>
    /// Асинхронный вспомогательный метод для обратного поиска
    /// </summary>
    /// <param name="strok"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <returns></returns>
    async private Task AsyncCalcualtion1Back(int strok, int i, int j)
    {
        D_backCopy[i] = await Add(D_backCopy[i], await Multipl(L_Matrix[strok, i], D_back[j], K_Path));
    }

    /// <summary>
    /// Сравнение
    /// </summary>
    /// <param name="D_back"></param>
    /// <param name="D_backCopy"></param>
    /// <returns></returns>
    private static async Task<Matrix<double>> Add(Matrix<double> D_back, Matrix<double> D_backCopy)
    {
        var matrixMIN = CreateMatrix.Dense<double>(1, D_back.ToArray().Length, 0);
        var Result_Matr = CreateMatrix.Dense(1, D_backCopy.ToArray().Length * 2, double.PositiveInfinity);
        for (int i = 0; i < D_back.ToArray().Length; i++)
        {
            Result_Matr[0, i] = D_back[0, i];
        }
        int iter = 0;
        for (int i = Result_Matr.ToArray().Length - 1; i >= Result_Matr.ToArray().Length - D_backCopy.ToArray().Length; i--)
        {
            Result_Matr[0, i] = D_backCopy[0, iter];
            iter++;
        }
        var array = Result_Matr.ToRowArrays();
        array[0] = array[0].OrderBy(x => x).ToArray();
        array[0] = array[0].Where(x => x != double.PositiveInfinity).ToArray();
        array[0] = array[0].Distinct().ToArray();
        for (int i = 0; i < matrixMIN.ToArray().Length; i++)
        {
            array[0] = array[0].Append(double.PositiveInfinity).ToArray();
        }
        Result_Matr = CreateMatrix.DenseOfRowArrays(array);
        array = null;
        iter = 0;
        for (int i = 0; i < matrixMIN.ToArray().Length; i++)
        {
            matrixMIN[0, i] = Result_Matr[0, i];
        }
        matrixMIN.ToRowArrays()[0] = matrixMIN.ToRowArrays()[0].OrderBy(x => x).ToArray();
        await Task.Delay(0);
        return matrixMIN;
    }
    /// <summary>
    /// Умножение
    /// </summary>
    /// <param name="matrix"></param>
    /// <param name="D"></param>
    /// <param name="K_Path"></param>
    /// <returns></returns>
    async private static Task<Matrix<double>> Multipl(Matrix<double> matrix, Matrix<double> D, int K_Path)
    {
        var matrixMIN = CreateMatrix.Dense<double>(1, K_Path, 0);
        int k = 0;
        Matrix<double> matrica = CreateMatrix.Dense<double>(1, K_Path * K_Path, 0);
        for (int i = 0; i < K_Path; i++)
        {
            for (int j = 0; j < K_Path; j++)
            {
                matrica[0, k] = Sum_INF(D[0, i], matrix[0, j]); //D[0, i] + matrix[0, j];
                k++;
            }
        }
        var array = matrica.ToRowArrays();
        array[0] = array[0].OrderBy(x => x).ToArray();
        matrica = CreateMatrix.DenseOfRowArrays(array);
        int iter = 0;

        while (iter != 3)
        {
            matrixMIN[0, iter] = matrica[0, iter];
            iter++;
        }
        await Task.Delay(0);
        return matrixMIN;
    }
    /// <summary>
    /// вспомогательная функция
    /// </summary>
    /// <param name="v1"></param>
    /// <param name="v2"></param>
    /// <returns></returns>
    private static double Sum_INF(dynamic v1, dynamic v2)
    {
        if (v1 == double.PositiveInfinity & v2 != double.PositiveInfinity)
            return v1;
        if (v1 != double.PositiveInfinity & v2 == double.PositiveInfinity)
            return v2;
        if (v1 == double.PositiveInfinity & v2 == double.PositiveInfinity)
            return double.PositiveInfinity;
        if (v1 != double.PositiveInfinity & v2 != double.PositiveInfinity)
            return v1 + v2;
        return 0;
    }

    /// <summary>
    /// Главный метод
    /// </summary>
    static int Main()
    {
        DoubleSearch doublesearch = new DoubleSearch();
        doublesearch.task = new Task(() => doublesearch.ALGORITHM());
        doublesearch.task.Start();
        doublesearch.task.Wait();
        Console.ReadKey();
        return 0;
    }
}

