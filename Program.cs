using System;
using Optimization;
using static System.Math;
using static Optimization.OptimizationMethods;

namespace Lab4
{
    class Program
    {

        static void Main(string[] args)
        {
            Function f;
            double answ = InputRequest("Выберите метод\n1 - Метод адаптивного поиска\n2 - Метод сверхбыстрого отжига\n");
            if(answ == 1)
            {
                answ = InputRequest("Выберите функцию\n" +
                "1 : 3 / x^2 + x\n" +
                "2 : -Sin(x) - Sin(3 * x) / 3\n" +
                "3 : (x - 1)^2)\n" +
                "4 : (x + 2.5) / (4 - x^2))\n" +
                "5 : 4 * x^3 - 8 * x^2 - 11 * x + 5\n" +
                "6 : Функция Химмельблау № 1\n" +
                "7 : Функция Химмельблау № 2\n" +
                "8 : Функция Вуда\n" +
                "9 : Функция Пауэлла");
                switch (answ)
                {
                    case 1:
                        f = new Function(f1, 1);
                        break;
                    case 2:
                        f = new Function(f2, 1);
                        break;
                    case 3:
                        f = new Function(f3, 1);
                        break;
                    case 4:
                        f = new Function(f4, 1);
                        break;
                    case 5:
                        f = new Function(f5, 1);
                        break;
                    case 6:
                        f = new Function(HimmelblauFunc1, 2);
                        break;
                    case 7:
                        f = new Function(HimmelblauFunc2, 2);
                        break;
                    case 8:
                        f = new Function(WoodFunc, 4);
                        break;
                    case 9:
                        f = new Function(PowellFunc, 4);
                        break;
                    default:
                        Console.WriteLine("Неизвестная функция, программа завершается");
                        return;
                }
                double[] startPos = new double[f.Rank];
                for (int i = 0; i < f.Rank; i++)
                {
                    startPos[i] = InputRequest($"Введите начальное значение x{i}");
                }
                double expansion = InputRequest("Введите коэффиициент расширения [1.1]");
                double compression = InputRequest("Введите коэффиициент сжатия [0.5]");
                int maxErrorCount = (int)InputRequest("Введите максимальное число ошибок на итерации [10]");
                double h = InputRequest("Введите начальную длину шага [0.5]");
                double R = InputRequest("Введите минимальную длину шага [0.0001]");
                int maxIterCount = (int)InputRequest("Введите максимальное число итераций [500]");

                double[] resX = MethodAdaptiveSearch(f, startPos, expansion, compression, maxErrorCount, h, R, maxIterCount);
                WriteABeautifullAnswer("По методу адаптивного поиска:", f[resX], resX);
            }
            else if(answ == 2)
            {
                answ = InputRequest("Выберите функцию\n" +
                "1 : Функция Растригина\n" +
                "2 : Функция Растригина овражная\n" +
                "3 : Функция Griewank\n" +
                "4 : Функция Катковника\n" +
                "5 : Функция Безымянная \n" +
                "6 : Функция Мультипликативная потенциальная\n");
                double S0, S1;
                switch (answ)
                {
                    case 1:
                        f = new Function(RastiginFunc, 2);
                        S0 = -16;
                        S1 = 16;
                        break;
                    case 2:
                        f = new Function(RastiginOvragFunc, 2);
                        S0 = -16;
                        S1 = 16;
                        break;
                    case 3:
                        f = new Function(GriewankFunc, 2);
                        S0 = -16;
                        S1 = 16;
                        break;
                    case 4:
                        f = new Function(KatkovnikFunc, 2);
                        S0 = -2.5;
                        S1 = 2.5;
                        break;
                    case 5:
                        f = new Function(SomeFunc9, 2);
                        S0 = 0;
                        S1 = 4;
                        break;
                    case 6:
                        f = new Function(MultiplicativePotentialFunc, 2);
                        S0 = 0;
                        S1 = 4;
                        break;
                    default:
                        Console.WriteLine("Неизвестная функция, программа завершается");
                        return;
                }
                double[] resX = MethodVeryFastAnnealin(f, S0, S1);
                WriteABeautifullAnswer("По методу сверхбыстрого отжига:", f[resX], resX);
            }
            else
            {
                return;
            }
            Console.ReadKey();
        }
        static double f1(params double[] x)
        {
            double value = 3 / Pow(x[0], 2) + x[0];
            return value;
        }
        static double f2(params double[] x)
        {
            double value = -Sin(x[0]) - Sin(3 * x[0]) / 3;
            return value;
        }
        static double f3(params double[] x)
        {
            double value = Pow(x[0] - 1, 2);
            return value;
        }
        static double f4(params double[] x)
        {
            double value = (x[0] + 2.5) / (4 - Pow(x[0], 2));
            return value;
        }
        static double f5(params double[] x)
        {
            double value = 4 * Pow(x[0], 3) - 8 * Pow(x[0], 2) - 11 * x[0] + 5;
            return value;
        }
        //Первая лаба многомерная оптимизация
        static double HimmelblauFunc1(params double[] x)//(5, 6) f = 0
        {
            double value = 4 * Pow((x[0] - 5), 2) + Pow(x[1] - 6, 2);
            return value;
        }
        static double HimmelblauFunc2(params double[] x)
        {
            double value = Pow(Pow(x[0], 2) + x[1] - 11, 2) + Pow(x[0] + Pow(x[1], 2) - 7, 2);//(3.58, - 1.85) f ~= 0.0011
            return value;
        }
        static double WoodFunc(params double[] x)//(-1.07, 1.16, -0.86, 0.76) f ~ 7.89
        {
            double value = 100 * Pow(x[1] - Pow(x[0], 2), 2) + Pow(1 - x[0], 2) + 90 *
                Pow(x[3] - Pow(x[2], 2), 2) + Pow(1 - x[2], 2) + 10.1 * (Pow(x[1] - 1, 2) +
                Pow(x[3] - 1, 2)) + 19.8 * (x[1] - 1) * (x[3] - 1);
            return value;
        }
        static double PowellFunc(params double[] x)//(0, 0, 0, 0) f = 0
        {
            double value = Pow(x[0] + 10 * x[1], 2) + 5 * Pow(x[2] - x[3], 2) + Pow(x[1] - 2 * x[2], 4) + 10 * Pow(x[0] - x[3], 4);
            return value;
        }

        static double SombreroFunc(params double[] x)//Сломано
        {
            double numerator = 1.0 - Pow(Sin(Sqrt(x[0] * x[0] + x[1] * x[1])), 2);
            double denumerator = 1.0 + 0.001 * (x[0] * x[0] + x[1] * x[1]);

            return numerator / denumerator;
        }
        static double DeJong2Func(params double[] x)//Сломано
        {
            return 100 - 100 / (100 * (x[0] * x[0] - x[1]) + Pow(1 - x[0], 2));
        }
        static double RastiginFunc(params double[] x)
        {
            return 0.1 * x[0] * x[0] + 0.1 * x[1] * x[1] - 4 * Cos(0.8 * x[0]) - 4 * Cos(0.8 * x[1]) + 8;
        }
        static double RastiginOvragFunc(params double[] x)
        {
            double Kx = 1.5;
            double Ky = 0.8;
            double a = Math.PI / 2;
            double A = x[0] * Cos(a) - x[1] * Sin(a);
            double B = x[0] * Sin(a) + x[1] * Cos(a);

            return Pow(0.1 * Kx * A, 2) + Pow(0.1 * Ky * B, 2) - 4 * Cos(0.8 * Kx * A) - 4 * Cos(0.8 * Ky * B) + 8;
        }
        static double GriewankFunc(params double[] x)
        {
            return 10 - 10 / (0.005 * (x[0] * x[0] + x[1] * x[1]) - Cos(x[0]) * Cos(x[1] / Sqrt(2)) + 2);
        }
        static double KatkovnikFunc(params double[] x)
        {
            double A = 0.8;
            return 0.5 * (Pow(x[0], 2) + Pow(x[1], 2)) * (2 * A + A * Cos(1.5 * x[0]) * Cos(Math.PI * x[1]) + A * Cos(Sqrt(5) * x[0]) * Cos(3.5 * x[1]));
        }
        static double SomeFunc9(params double[] x)
        {
            return 0.5 * (x[0] * x[0] + x[0] * x[1] + x[1] * x[1]) * (1 + 0.5 * Cos(1.5 * x[0]) * Cos(3.2 * x[0] * x[1]) * Cos(3.14 * x[1]) + 0.5 * Cos(2.2 * x[0]) * Cos(4.8 * x[0] * x[1]) * Cos(3.5 * x[1]));
        }
        static double MultiplicativePotentialFunc(params double[] x)
        {
            return - zforMultiplicativePotentialFunc(x[1]) * zforMultiplicativePotentialFunc(x[0]);
        }
        static double zforMultiplicativePotentialFunc(double x)
        {
            return -(1.0 / (Pow(x - 1, 2) + 0.2)) - (1.0 / (2 * Pow(x - 2, 2) + 0.15)) - (1.0 / (3 * Pow(x - 3, 2) + 0.3));
        }


    }
}






//int k = 0;
//double[] startPos = new double[2];
//double S0 = 0;
//double S1 = 4;
//            for (int i = 0; i< 2; i++)
//            {
//                startPos[i] = random.NextDouble() * (S1 - S0) + S0;
//            }
//            Function f = new Function(MultiplicativePotentialFunc, 2);
//double Energy = f[startPos];
//double[] curX = (double[])startPos.Clone();
//double[] tryX;
//double H = 0;
//double deltaE = 0;
//double[] GlobalMin = (double[])startPos.Clone();
//double curT = 20;
//            while (k< 1000000)
//            {
//                if (Energy<f[GlobalMin])
//                {
//                    GlobalMin = (double[]) curX.Clone();

//                }
//                double curRand = 0;
//                do
//                {
//                    tryX = nextX(curX, S0, S1, k);
//deltaE = f[tryX] - E;
//                    curT = DamnT(k, 2);
//H = DamnH(deltaE, curT);
//curRand = random.NextDouble();
//                    double alphaMinusH = curRand - H;
//                } while (H <= curRand);
//                curX = (double[]) tryX.Clone();
//Energy = f[curX];
//                k++;
//                //WriteABeautifullAnswer($"Текущая итерация = {k} A Энергия = {Energy}", f[GlobalMin], GlobalMin);
//            }
//            WriteABeautifullAnswer($"Последняя итерация = {k} A Энергия = {f[GlobalMin]}", f[GlobalMin], GlobalMin);