using System;
using System.Linq;
using static System.Math;

namespace Optimization
{
    static public class OptimizationMethods
    {
        public static Random random = new Random();
        /// <summary>
        /// Метод Дэвиса-Свенна-Кэмпи для поиск отрезка минимума [a,b] 
        /// </summary>
        /// <param name="f"> Минимизируемая функция </param>
        /// <param name="h"> Длина шага </param>
        /// <param name="a"> Начало отрезка минимума функции </param>
        /// <param name="b"> Конец отрезка минимума функции </param>
        /// <param name="startPos"> Точка начала поиска </param>
        /// <param name="searchDirection"> Направление поиска для функции двух и более переменнных</param>
        /// <param name="xNum"> Номер аргумента минимизации для функции двух и более переменнных</param>
        /// <param name="x0"> Начальные значения минимизации для функции двух и более переменнных</param>
        static public void MethodDSK(Function f, double h, out double a, out double b, double startPos, int xNum = 0, double searchDirection = 0, params double[] x0)
        {
            double curX = startPos;
            int rank = f.Rank;
            double[] fVal = new double[3];//значения функций 
            if (rank == 1)//Функция минимизации - одномерная функция
            {
                fVal[0] = f[curX - h];
                fVal[1] = f[curX];
                fVal[2] = f[curX + h];
            }
            else//Функция минимизации - многомерная функция
            {
                double[] currentXs = (double[])x0.Clone();//Массив текущих значений иксов 
                currentXs[xNum] = x0[xNum] + (curX - h) * searchDirection;
                fVal[0] = f[currentXs];
                currentXs[xNum] = x0[xNum] + curX * searchDirection;
                fVal[1] = f[currentXs];
                currentXs[xNum] = x0[xNum] + (curX + h) * searchDirection;
                fVal[2] = f[currentXs];
            }

            if (fVal[0] >= fVal[1] && fVal[2] >= fVal[1])//Если значения функции слева и справа больше значения функции в начальной точке
            {
                a = curX - h;
                b = curX + h;
                return;
            }
            else if (fVal[2] > fVal[1])//Если значения функции убывают слева направо, то продолжить поиск
            {
                h = -h;
                curX += h;
            }
            else
            {
                curX += h;
            }
            int k = 0;
            Array.Resize(ref fVal, 2);
            do
            {
                if (rank == 1)
                {
                    fVal[0] = f[curX - h];
                    fVal[1] = f[curX];
                }
                else
                {
                    double[] currentXs = (double[])x0.Clone();
                    currentXs[xNum] = x0[xNum] + (curX - h) * searchDirection;
                    fVal[0] = f[currentXs];
                    currentXs[xNum] = x0[xNum] + curX * searchDirection;
                    fVal[1] = f[currentXs];
                }
                curX += h;
                k++;
            } while (fVal[1] < fVal[0] && k < 100000);
            curX -= h;
            if (h > 0)
            {
                a = curX - 2 * h;
                b = curX;
            }
            else
            {
                a = curX;
                b = curX - 2 * h;
            }
        }
        /// <summary>
        /// /Метот параболической аппроксимации Пауэла 
        /// </summary>
        /// <param name="f"> Минимизируемая функция </param>
        /// <param name="a"> Начало отрезка минимума функции </param>
        /// <param name="b"> Конец отрезка минимума функции </param>
        /// <param name="eps"> Точность минимизации </param>
        /// <param name="searchDirection"> Направление поиска для функции двух и более переменнных</param>
        /// <param name="xNum"> Номер аргумента минимизации для функции двух и более переменнных</param>
        /// <param name="x0"> Начальные значения минимизации для функции двух и более переменнных</param>
        /// <returns>Точку минимума функции</returns>
        static public double MethodParabolicApproximation(Function f, double a, double b, double eps = 0.01, double searchDirection = 0, int xNum = 0, params double[] x0)
        {
            int rank = f.Rank;
            double[] x = (double[])x0.Clone();
            double[] z = new double[4];
            z[0] = a;
            z[1] = a / 2 + b / 2;
            z[2] = b;
            double zAvr;
            double zMin;
            int k = 0;
            do
            {
                double[] fVal = new double[4];
                if (rank == 1)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        fVal[i] = f[z[i]];
                    }
                }
                else
                {
                    double[] currentXs = (double[])x0.Clone();
                    for (int i = 0; i < 3; i++)
                    {
                        currentXs[xNum] = x[xNum] + z[i] * searchDirection;
                        fVal[i] = f[currentXs];
                    }
                }
                z[3] = (Pow(z[1], 2) - Pow(z[2], 2)) * fVal[0] + (Pow(z[2], 2) - Pow(z[0], 2)) * fVal[1] + (Pow(z[0], 2) - Pow(z[1], 2)) * fVal[2];
                z[3] /= 2;
                z[3] /= (z[1] - z[2]) * fVal[0] + (z[2] - z[0]) * fVal[1] + (z[0] - z[1]) * fVal[2];
                if (double.IsNaN(z[3]) || double.IsInfinity(z[3]))
                {
                    z[3] = Max(Max(z[1], z[2]), z[0]);
                }    


                if (rank == 1) fVal[3] = f[z[3]];
                else
                {
                    double[] currentXs = (double[])x0.Clone();
                    currentXs[xNum] = x[xNum] + z[3] * searchDirection;
                    fVal[3] = f[currentXs];
                }
                zAvr = z[3];
                if (fVal[0] <= fVal[1] && fVal[0] <= fVal[2])
                {
                    zMin = z[0];
                }
                else if (fVal[1] <= fVal[0] && fVal[1] <= fVal[2])
                {
                    zMin = z[1];
                }
                else
                {
                    zMin = z[2];
                }
                double fMax = Max(Max(fVal[0], fVal[1]), Max(fVal[2], fVal[3]));
                double[] tmpZ = new double[4];
                bool maxZisFound = false;
                for (int i = 0; i < 4; i++)
                {
                    if ((fVal[i]) != fMax || maxZisFound)
                    {
                        tmpZ[i] = z[i];
                    }
                    else
                    {
                        tmpZ[i] = double.PositiveInfinity;
                        maxZisFound = true;
                    }
                }
                Array.Copy(tmpZ, z, 4);
                Array.Sort(z);
                k++;
            }
            while ((Abs(zAvr - zMin) > eps && k < 100));

            return zAvr;
        }
        /// <summary>
        /// Метод Дэвидона-Флетчера-Пауэлла 
        /// </summary>
        /// <param name="f"> Минимизируемая функция </param>
        /// <param name="startPos"> Массив начальных значений функции </param>
        /// <param name="eps"> Точность минимизации </param>
        /// <returns></returns>
        static public double[] MethodDFP(Function f, double[] startPos, double eps)//Дэвидона - Флетчра - Пауэла
        {
            int k = 0;
            int rank = f.Rank;
            double[] gradArr;
            double[] curX = (double[])startPos.Clone();
            double[,] hessian = new double[rank, rank];//Матрица Гессе
            for (int i = 0; i < rank; i++)
            {
                hessian[i, i] = 1;
            }
            double lymbda; //, lymbdaStart, lymbdaEnd;
            double[,] helpVectS = new double[rank, rank];//Вспомогательный вектор
            double[,] helpVectY = new double[rank, rank];
            while (k < 500)
            {
                gradArr = Gradient(f, curX);
                if (EuclideanNorm(gradArr) < eps)
                {
                    return curX;
                }
                if (k != 0)
                {
                    double[] prevStartPos = (double[])curX.Clone();
                    for (int i = 0; i < rank; i++)
                    {
                        prevStartPos[i] -= helpVectS[i, 0];
                    }
                    double[] prevGradArr = Gradient(f, prevStartPos);
                    for (int i = 0; i < rank; i++)
                    {
                        helpVectY[i, 0] = gradArr[i] - prevGradArr[i];
                    }
                    //Числитель 1
                    double[,] numerator1 = MatrixMul(MatrixMul(MatrixMul(hessian, helpVectY), MatrixTranspose(helpVectY)), hessian);
                    //Знаменатель 1
                    double denominator1 = MatrixMul(MatrixMul(MatrixTranspose(helpVectY), hessian), helpVectY)[0, 0];
                    double[,] buffer1 = MatrixOnConst(numerator1, 1.0 / denominator1);
                    //Числитель 2
                    double[,] numerator2 = MatrixMul(helpVectS, MatrixTranspose(helpVectS));
                    //Знаменатель 2
                    double denominator2 = MatrixMul(MatrixTranspose(helpVectY), helpVectS)[0, 0];
                    double[,] buffer2 = MatrixOnConst(numerator2, 1.0 / denominator2);

                    for (int i = 0; i < rank; i++)
                    {
                        for (int j = 0; j < rank; j++)
                        {
                            hessian[i, j] -= buffer1[i, j];
                            hessian[i, j] += buffer2[i, j];
                        }
                    }
                }

                double[] searchDirection = new double[rank];//Квазиньютоновское направление поиска
                for (int i = 0; i < rank; i++)
                {
                    for (int j = 0; j < rank; j++)
                    {
                        searchDirection[i] += -hessian[i, j] * gradArr[i];
                    }
                }
                for (int i = 0; i < rank; i++)
                {
                    if (searchDirection[i] != 0)
                    {
                        double[] z;
                        lymbda = 2;
                        do
                        {
                            lymbda /= 2;
                            z = (double[])curX.Clone();
                            z[i] += lymbda * searchDirection[i];
                        } while (Abs(f[z] - f[curX]) > eps && (f[z] > f[curX]));
                        helpVectS[i, 0] = searchDirection[i] * lymbda;
                        curX[i] += helpVectS[i, 0];
                    }
                }
                k++;
            }
            Console.WriteLine(string.Format("Превышено число итераций, ответ может быть неверен"));
            return curX;
        }
        static public double[] MethodNR(Function f, double[] startPos, double eps)//Метод Ньютона-Рафсона
        {
            int k = 0;
            int rank = f.Rank;
            double[] gradVector;
            double[] curX = (double[])startPos.Clone();
            double[,] hessian;//Матрица Гессе
            double lymbda;
            while (k < 1000)
            {
                gradVector = Gradient(f, curX);
                if (EuclideanNorm(gradVector) <= eps)
                {
                    return curX;
                }

                hessian = FillHessian(f, curX);
                double[,] inverceHessian = MatrixInverse(hessian);
                for (int i = 0; i < rank; i++)
                {
                    for (int j = 0; j < rank; j++)
                    {
                        inverceHessian[i, j] = -inverceHessian[i, j];
                    }
                }

                double[,] searchDirection;//Направление поиска

                double[,] gradForMatixMul = new double[rank, rank];
                for (int i = 0; i < rank; i++)
                {
                    gradForMatixMul[i, 0] = gradVector[i];
                }

                searchDirection = MatrixMul(inverceHessian, gradForMatixMul);

                for (int i = 0; i < rank; i++)
                {
                    if (searchDirection[i, 0] != 0)
                    {
                        double[] z;
                        lymbda = 2;
                        do
                        {
                            lymbda /= 2;
                            z = (double[])curX.Clone();
                            z[i] += lymbda * searchDirection[i, 0];
                        } while (Abs(f[z] - f[curX]) > eps && (f[z] > f[curX]));
                        curX[i] += lymbda * searchDirection[i, 0];
                    }
                }
                k++;
            }
            Console.WriteLine(string.Format("Превышено число итераций, ответ может быть неверен"));
            return curX;
        }
        /// <summary>
        /// Метод Гаусса - Зейделя
        /// </summary>
        /// <param name="f"> Минимизируемая функция </param>
        /// <param name="curPos"> Массив начальных значений функции </param>
        /// <param name="h"> Массив длины шагов для каджого параметра функции </param>
        /// <param name="eps"> Точность минимизации </param>
        /// <returns></returns>
        static public double[] MethodGS(Function f, double[] startPos, double[] h, double eps)
        {
            int rank = f.Rank;
            int k = 0;
            double lymbdaStart, lymbda , lymbdaEnd;
            double z = 0.5;
            double[] curPos = (double[])startPos.Clone();
            do
            {
                double[] prevX = (double[])curPos.Clone(); 
                for (int i = 0; i < rank; i++)
                {
                    MethodDSK(f, h[i], out lymbdaStart, out lymbdaEnd, 0, i, h[i], curPos);
                    lymbda = MethodParabolicApproximation(f, lymbdaStart, lymbdaEnd, eps, h[i], i, curPos);
                    lymbda= Round(lymbda);
                    curPos[i] += lymbda * h[i];
                }
                if(Enumerable.SequenceEqual(curPos, prevX))
                {
                   if(EuclideanNorm(h) <= eps) 
                        return curPos;
                   for (int j = 0; j < rank; j++) 
                        h[j] *= z;

                }
                k++;
            } while (k < 100 && EuclideanNorm(h) > eps);
            if(k == 100)
            {
                Console.WriteLine("Превышено число итераций, ответ может быть неверным");
            }
            return curPos;
        }
        static public double[] MethodPenaltyFunctions(Function f,BarrierFunction bf, PenaltyFunction pf,
            double[] startPos, double startPenalty, double mulFactor = 4, double eps = 0.01)
        {
            int k = 0;
            double penalty = startPenalty;
            double[] curX = (double[])startPos.Clone();
            while ( k < 100)
            {
                Function omg = new Function(pf, bf, f, penalty);
                curX = MethodNR(omg, curX, eps);
                if(penalty * (bf[curX] + pf[curX]) < eps)
                {
                    return curX;
                }
                penalty *= mulFactor;
                k++;
            }
            return curX;
        }
        static public double[] MethodBarrierFunctions(Function f, BarrierFunction bf,double[] startPos, double startPenalty = 1, double mulFactor = 0.1,  double eps = 0.01)
        {
            int k = 0;
            double penalty = startPenalty;
            double[] curX = (double[])startPos.Clone();
            while (k < 100)
            {
                Function omg = new Function(bf, f, penalty);
                curX = MethodNR(omg, curX, eps);
                if (Abs(penalty / bf[curX]) < eps)
                {
                    return curX;
                }
                penalty *= mulFactor;
                k++;
            }
            return curX;
        }
        static public double[] MethodAdaptiveSearch(Function f, double[] startPos, double expansion = 1.1, double compression = 0.5, int maxErrorCount = 10, double h = 0.5, double R = 0.0001, int maxIterCount = 500)
        {
            int Rank = f.Rank;
            int k = 0;//Номер шага
            int j = 1;
            double[] randVect = new double[Rank];
            Random random = new Random();
            double[] curX = (double[])startPos.Clone();
            bool wrongSeachDirection = false;
            while (true)
            {
                for (int i = 0; i < f.Rank; i++)
                {
                    randVect[i] = random.NextDouble() * 2 - 1;
                }
                double[] nextX = new double[Rank];
                for (int i = 0; i < f.Rank; i++)
                {
                    nextX[i] = curX[i] + h * randVect[i] / EuclideanNorm(randVect);
                }
                if (f[nextX] < f[curX] && !double.IsInfinity(f[curX]))
                {
                    double[] z = new double[Rank];
                    for (int i = 0; i < f.Rank; i++)
                    {
                        z[i] = curX[i] + expansion * (nextX[i] - curX[i]);
                    }
                    if (f[z] < f[curX])
                    {
                        curX = (double[])z.Clone();
                        h *= expansion;
                        k++;
                        if (k >= maxIterCount)
                        {
                            return curX;
                        }
                        else
                        {
                            j = 1;
                        }
                    }
                    else//Шаг 5
                    {
                        wrongSeachDirection = true;
                    }
                }
                else//Шаг 5
                {
                    wrongSeachDirection = true;
                }
                if (wrongSeachDirection)
                {
                    if (j < maxErrorCount)
                    {
                        j++;
                    }
                    else
                    {
                        if (h <= R)
                        {
                            return curX;
                        }
                        else
                        {
                            h *= compression;
                            j = 1;
                        }
                    }
                    wrongSeachDirection = false;
                }
            }
        }
        static public double[] MethodVeryFastAnnealin(Function f, double S0, double S1)
        {
            int k = 0;
            double[] startPos = new double[2];
            for (int i = 0; i < 2; i++)
            {
                startPos[i] = random.NextDouble() * (S1 - S0) + S0;
            }
            double Energy = f[startPos];
            double[] curX = (double[])startPos.Clone();
            double[] tryX;
            double H = 0;
            double deltaE = 0;
            double[] GlobalMin = (double[])startPos.Clone();
            double curT = 20;
            while (k < 1000000)
            {
                if (Energy < f[GlobalMin])
                {
                    GlobalMin = (double[])curX.Clone();

                }
                double curRand = 0;
                do
                {
                    tryX = nextX(curX, S0, S1, k);
                    deltaE = f[tryX] - E;
                    curT = DamnT(k, 2);
                    H = DamnH(deltaE, curT);
                    curRand = random.NextDouble();
                    double alphaMinusH = curRand - H;
                } while (H <= curRand);
                curX = (double[])tryX.Clone();
                Energy = f[curX];
                k++;
            }
            return GlobalMin;
        }
        static double[] nextX(double[] x0, double S0, double S1, int k)
        {
            double[] curX = new double[x0.Length];
            for (int i = 0; i < x0.Length; i++)
            {
                do
                {
                    double curT = DamnT(k, 2);

                    curX[i] = x0[i] + DamnZ(curT) * (S1 - S0);
                } while (curX[i] > S1 || curX[i] < S0);
            }
            return curX;
        }
        static double DamnZ(double T)
        {
            double rand = random.NextDouble();
            double sgn = Sign(rand - 0.5);
            double power = Abs(2 * rand - 1);

            return sgn * T * (Pow(1 + (1 / T), power) - 1);
        }
        static double DamnT(int k, int D)
        {
            double T0 = 5;
            double ci = 0.5;
            double kPow = Pow(k, 1.0 / D);
            double answ = T0 * Exp(-ci * kPow);
            return answ;
        }
        static double DamnH(double deltaE, double T)
        {
            return Exp(-deltaE / T);
        }

        /// <summary>
        /// Функция запроса ввода числа 
        /// </summary>
        /// <param name="requestMessage"></param>
        /// <returns></returns> 
        public static double InputRequest(string requestMessage = "Введите число")
        {
            double value;
            Console.WriteLine(requestMessage);
            while(!double.TryParse(Console.ReadLine(), out value))
            {
                Console.WriteLine("Введённое значение некоректно, повторите попытку ввода");
            }
            Console.WriteLine("Вы ввели значение равное {0} ", value);
            return value;
        }
        public static double[,] MatrixMul(double[,] matrix1, double[,] matrix2)//Произведение двух матриц
        {
            if (matrix1.GetLength(0) != matrix2.GetLength(1))
            {
                Console.WriteLine("Неверные размерности матриц");
                return new double[matrix1.GetLength(0), matrix1.GetLength(1)];
            }
            double[,] resMatrix = new double[matrix1.GetLength(1), matrix2.GetLength(0)];
            for (int i = 0; i < matrix1.GetLength(0); i++)
            {
                for (int j = 0; j < matrix2.GetLength(1); j++)
                {
                    for (int k = 0; k < matrix2.GetLength(0); k++)
                    {
                        resMatrix[i, j] += matrix1[i, k] * matrix2[k, j];
                    }
                }
            }
            return resMatrix;
        }
        public static double[,] MatrixTranspose(double[,] matrix)//Транспонирование матрицы
        {
            double[,] buffer = (double[,])matrix.Clone();
            for (int i = 0; i < buffer.GetLength(0); i++)
            {
                for (int j = i; j < buffer.GetLength(1); j++)
                {
                    double buff = buffer[i, j];
                    buffer[i, j] = matrix[j, i];
                    buffer[j, i] = buff;
                }
            }
            return buffer;
        }
        public static double[,] MatrixInverse(double[,] matrix)//Нахождение обратно мартицы
        {
            int rank = matrix.GetLength(0);
            double[,] reversedMatrix = (double[,])matrix.Clone();
            double det = MatrixDeterminant(reversedMatrix);
            double[,] minorMarix = new double[rank, rank];
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < rank; j++)
                {
                    double[,] minor = new double[rank - 1, rank - 1];
                    for (int a = 0, a1 = 0; a < rank; a++)
                    {
                        int b1 = 0;
                        for (int b = 0; b < rank; b++)
                        {
                            if (a != i && b != j) {
                                minor[a1, b1] = reversedMatrix[a, b];
                                b1++;
                            }
                        }
                        if (a != i)
                            a1++;
                    }
                    minorMarix[i, j] = Pow(-1.0, i + j) *MatrixDeterminant(minor);
                }
            }
            reversedMatrix = MatrixTranspose(minorMarix);
            reversedMatrix = MatrixOnConst(reversedMatrix, 1.0 / det);

            return reversedMatrix;
        }
        public static double MatrixDeterminant(double[,] matrix)//Нахождение определителя матрицы 
        {

            if(matrix == null)
            {
                throw new ArgumentNullException("MatrixDeterminant(double[,] matrix) = null");
            }
            double det;
            if(matrix.GetLength(0) == 1)
            {
                det = matrix[0, 0];
            }
            else if (matrix.GetLength(0) == 2)
            {
                det = (matrix[0, 0] * matrix[1, 1]) - (matrix[0, 1] * matrix[1, 0]);
            }
            else if(matrix.GetLength(0) == 3)
            {
                det =
                matrix[0, 0] * matrix[1, 1] * matrix[2, 2] +
                matrix[0, 1] * matrix[1, 2] * matrix[2, 0] +
                matrix[0, 2] * matrix[1, 0] * matrix[2, 1] -
                matrix[0, 2] * matrix[1, 1] * matrix[2, 0] -
                matrix[0, 0] * matrix[1, 2] * matrix[2, 1] -
                matrix[0, 1] * matrix[1, 0] * matrix[2, 2];

            }
            else if(matrix.GetLength(0) == 4)
            {
                det = 
                matrix[0, 3] * matrix[1, 2] * matrix[2, 1] * matrix[3, 0] - matrix[0, 2] * matrix[1, 3] * matrix[2, 1] * matrix[3, 0] -
                matrix[0, 3] * matrix[1, 1] * matrix[2, 2] * matrix[3, 0] + matrix[0, 1] * matrix[1, 3] * matrix[2, 2] * matrix[3, 0] +
                matrix[0, 2] * matrix[1, 1] * matrix[2, 3] * matrix[3, 0] - matrix[0, 1] * matrix[1, 2] * matrix[2, 3] * matrix[3, 0] -
                matrix[0, 3] * matrix[1, 2] * matrix[2, 0] * matrix[3, 1] + matrix[0, 2] * matrix[1, 3] * matrix[2, 0] * matrix[3, 1] +
                matrix[0, 3] * matrix[1, 0] * matrix[2, 2] * matrix[3, 1] - matrix[0, 0] * matrix[1, 3] * matrix[2, 2] * matrix[3, 1] -
                matrix[0, 2] * matrix[1, 0] * matrix[2, 3] * matrix[3, 1] + matrix[0, 0] * matrix[1, 2] * matrix[2, 3] * matrix[3, 1] +
                matrix[0, 3] * matrix[1, 1] * matrix[2, 0] * matrix[3, 2] - matrix[0, 1] * matrix[1, 3] * matrix[2, 0] * matrix[3, 2] -
                matrix[0, 3] * matrix[1, 0] * matrix[2, 1] * matrix[3, 2] + matrix[0, 0] * matrix[1, 3] * matrix[2, 1] * matrix[3, 2] +
                matrix[0, 1] * matrix[1, 0] * matrix[2, 3] * matrix[3, 2] - matrix[0, 0] * matrix[1, 1] * matrix[2, 3] * matrix[3, 2] -
                matrix[0, 2] * matrix[1, 1] * matrix[2, 0] * matrix[3, 3] + matrix[0, 1] * matrix[1, 2] * matrix[2, 0] * matrix[3, 3] +
                matrix[0, 2] * matrix[1, 0] * matrix[2, 1] * matrix[3, 3] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] * matrix[3, 3] -
                matrix[0, 1] * matrix[1, 0] * matrix[2, 2] * matrix[3, 3] + matrix[0, 0] * matrix[1, 1] * matrix[2, 2] * matrix[3, 3];
            }
            else
            {
                det = 0;
            }
            return det;
        }
        public static double[,] MatrixOnConst(double[,] matrix, double constValue)//Нахождение произведения матрицы на константу
        {
            double[,] buffer = new double[matrix.GetLength(0), matrix.GetLength(1)]; 
            for (int i = 0; i < buffer.GetLength(0); i++)
            {
                for (int j = 0; j < buffer.GetLength(1); j++)
                {
                    buffer[i, j] = matrix[i, j] * constValue;
                }
            }
            return buffer;
        }
        public static void WriteABeautifullAnswer(string msg,double functionValue, params double[] values)
        {
            Console.WriteLine(msg);
            string text = "";
            string text2 = "|";
            int count = 0;
            foreach(double val in values)
            {
                count += val.ToString("F4").Length;
            }
            count += functionValue.ToString("F4").Length;
            int spacesPerX = count / (values.Length + 1);
            for (int i = 0;i < count + 3 * values.Length + 6; i++)
            {
                text += "-";
            }
            
            Console.WriteLine(text);
            for (int i = 0; i < ((count + 3 * values.Length + 6) / 2)  - (1 + "Минимум функции найден".Length / 2); i++)
            {
                text2 += " ";
            }
            text2 += "Минимум функции найден";
            for (int i = 0; i < ((count + 3 * values.Length + 6) / 2) - (1 + "Минимум функции найден".Length / 2); i++)
            {
                text2 += " ";
            }
            text2 += "|";
            Console.WriteLine(text2);
            Console.WriteLine(text);
            string text3 = "| ";
            for(int i = 0; i < values.Length; i++)
            {
                for(int j = 0; j < (spacesPerX / 2) - 1; j++)
                {
                    text3 += " ";
                }
                text3 += string.Format("x{0}", i);
                for (int j = 0; j < (spacesPerX / 2); j++)
                {
                    text3 += " ";
                }
                text3 += "|";
            }
            for (int j = 0; j < (spacesPerX / 2) - 1; j++)
            {
                text3 += " ";
            }
            text3 += " f";
            for (int j = 0; j < (spacesPerX / 2); j++)
            {
                text3 += " ";
            }
            text3 += " |";

            Console.WriteLine(text3);
            Console.WriteLine(text);
            string text4 = "|";
            for (int i = 0; i < values.Length; i++)
            {
                for (int j = 0; j < (spacesPerX / 2) - values[i].ToString("F4").Length / 2; j++)
                {
                    text4 += " ";
                }
                text4 += values[i].ToString("F4");
                for (int j = 0; j < (spacesPerX / 2) - values[i].ToString("F4").Length / 2; j++)
                {
                    text4 += " ";
                }
                text4 += " | ";
            }
            for (int j = 0; j < (spacesPerX / 2) - functionValue.ToString("F4").Length / 2; j++)
            {
                text4 += " ";
            }
            text4 += functionValue.ToString("F4");
            for (int j = 0; j < (spacesPerX / 2) - functionValue.ToString("F4").Length / 2; j++)
            {
                text4 += " ";
            }
            text4 += " |";
            Console.WriteLine(text4);
            Console.WriteLine(text);
        }
        public static double[] Gradient(Function f, double[] x)
        {
            double delta = 1e-8;
            int rank = f.Rank;
            double[] gradArr = new double[rank];
            double[] curX = (double[])x.Clone();
            for (int i = 0; i < rank; i++)
            {
                curX[i] = x[i] + delta;
                gradArr[i] = (f[curX] - f[x]) / delta;
                curX[i] = x[i];
            }
            return gradArr;
        }
        public static double EuclideanNorm(params double[] vals)
        {
            double norm = 0;
            foreach (double val in vals) norm += Pow(val, 2);
            return Sqrt(norm);
        }
        public static double[,] FillHessian(Function f, double[] x) 
        {
            int rank = f.Rank;
            double[,] hessian = new double[rank, rank];
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < rank; j++)
                {
                    hessian[i, j] = SecondPatrialDerivative(f, x, i, j);
                }
            }
            return hessian;
        }
        public static double SecondPatrialDerivative(Function f, double[] x, int dx1, int dx2) 
        {
            double result;
            double delta = 1e-4;
            double[] curX = (double[])x.Clone();
            if (dx1 == dx2) 
                curX[dx1] = x[dx1] + 2 * delta;
            else {
                curX[dx1] = x[dx1] + delta;
                curX[dx2] = x[dx2] + delta;
            }
            double result1 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1];
            else
            {
                curX[dx1] = x[dx1] + delta;
                curX[dx2] = x[dx2] - delta;
            }
            double result2 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1];
            else
            {
                curX[dx1] = x[dx1] - delta;
                curX[dx2] = x[dx2] + delta;
            }
            double result3 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1] - 2 * delta;
            else
            {    
                curX[dx2] = x[dx2] - delta;   
            }
            double result4 = f[curX];
            result = (result1 - result2 - result3 + result4) / (4 * Pow(delta, 2));
            return result;
        }
    }
    public class Function
    {
        public delegate double function(params double[] values);
        readonly int type = 0;
        private readonly double penalty;
        private readonly PenaltyFunction pf;
        private readonly BarrierFunction bf;
        private readonly Function ff;
        public Function(function f, int rank)
        {
            this.f = f;
            Rank = rank;
        }
        public Function(PenaltyFunction pf, BarrierFunction bf, Function ff, double penalty)
        {
            type = 1;
            this.pf = pf;
            this.bf = bf;
            this.ff = ff;
            Rank = ff.Rank;
            this.penalty = penalty;
        }
        public Function(BarrierFunction bf, Function ff, double penalty)
        {
            type = 2;
            this.bf = bf;
            this.ff = ff;
            Rank = ff.Rank;
            this.penalty = penalty;
        }
        public int Rank { get; set; }//Размерность функции
        protected function f;
        virtual public double this[params double[] values]
        {
            get 
            {   if(type == 0)
                {
                    return f(values);
                }
                else if(type == 1)
                {
                    return penalty * (pf[values] + bf[values]) + ff[values];
                }
                else
                {
                    return ff[values] + ( - penalty / bf[values]);
                }

                
            }
        }

    }
    public class BarrierFunction : Function
    {
        public BarrierFunction(function f, int rank) : base(f, rank)
        {
            this.f = f;
            Rank = rank;
        }
        override public double this[params double[] values]
        {
            get 
            {
                return f(values);
            }
        }

    }
    public class PenaltyFunction: Function
    {
        public PenaltyFunction(function f, int rank) : base(f, rank)
        {
            this.f = f;
            Rank = rank;
        }
        override public double this[params double[] values]
        {
            get
            {
                return Pow(f(values), 2);
            }
        }

    }
}