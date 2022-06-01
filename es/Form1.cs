using Extreme.Mathematics;
using ILNumerics;
using ILNumerics.Drawing;
using ILNumerics.Drawing.Plotting;
using ILNumerics.Toolboxes;
using KdTreeDemo;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using static ILNumerics.Globals;
using static ILNumerics.ILMath;
using pts = KdTreeDemo.Points;
namespace es
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        /// <summary>
        /// 返回网格点
        /// </summary>
        /// <param name="scatteredValues">采样值</param>
        /// <param name="scatteredPositions">X,Y采样点</param>
        /// <param name="start">网格最小值</param>
        /// <param name="end">网格最大值</param>
        /// <param name="step">网格精度</param>
        /// <returns></returns>
        private Array<double> Kriging(Array<Double> scatteredValues, Array<double> scatteredPositions, double start, double end, float step)
        {

            using (Scope.Enter(ArrayStyles.ILNumericsV4))
            {
                Array<double> outPut = new double[3, (int)(Math.Pow(((Math.Abs(end) + Math.Abs(end)) / step), 2) + 1)];
                /*
                 * 创建网格
                 */
                Array<double> Y = 1, X = meshgrid(arange(start, step, end), arange(start, step, end), Y);
                Array<double> interpPositions = X[full].T.Concat(Y[full].T, 0);
                /*
                 * 克里金插值模型
                 */
                Array<double> err = 1;
                Array<double> interpValues = Interpolation.kriging(scatteredValues, scatteredPositions, interpPositions, error: err);
                // let's plot it! 
                var scene = new Scene();
                var pc = scene.Add(new PlotCube(twoDMode: false));
                /*
                 * 绘制点
                 */
                Array<float> interpValues4Plot = tosingle(scatteredPositions.Concat(scatteredValues, 0));
                pc.Add(new ILNumerics.Drawing.Points()
                {
                    Positions = interpValues4Plot,
                    Color = Color.LightPink,
                    Size = 7
                });
                /*
                 * 将这些值恢复到合适的形状以进行绘图
                 */
                interpValues.a = reshape(interpValues, X.shape);
                err.a = reshape(err, X.shape);
                outPut = X[":"].T.Concat(Y[":"].T, 0).Concat(interpValues[":"].T, 0);
                /*
                 * 绘制曲面
                 */
                pc.Add(new Surface(interpValues, X, Y, colormap: Colormaps.Jet, C: err)
                {
                    Fill = { Visible = false },
                    Wireframe = { Width = 2, Color = null }, // 'Color' controls the solid color. 'null' activates individual color mapped colors. 
                    Children = {
                    new Colorbar() {
                        // give the colorbar a title
                        new ILNumerics.Drawing.Label(text: @"Error, (\vartheta)") {
                            Anchor = new PointF(-.4f,1.1f)
                        }
                    }
                }
                });
                panel1.Scene = scene;
                // don't forget to prepare the scene for rendering + trigger a redraw!
                panel1.Configure();
                panel1.Refresh();
                return outPut;
            }
        }
        /// <summary>
        /// 返回目标点z值
        /// </summary>
        /// <param name="scatteredValues">采样点值</param>
        /// <param name="scatteredPositions">采样点</param>
        /// <param name="target">目标点</param>
        /// <param name="k">邻近K点</param>
        /// <returns></returns>
        private double Kriging(Array<Double> scatteredValues, Array<double> scatteredPositions, PointF target, int k)
        {
            double out_value = 0;
            using (Scope.Enter(ArrayStyles.ILNumericsV4))
            {
                /*
                 * 确定精度大小
                 */
                int target_x = target.X.ToString().Length - target.X.ToString().IndexOf(".");
                int target_y = target.Y.ToString().Length - target.Y.ToString().IndexOf(".");
                int retain = Math.Max(target_x, target_y) - 1;
                decimal step = DecimalMath.Pow((decimal)0.1, retain);
                /*
                 * 邻近K点
                 */
                List<pts> kPoints = new List<pts>();
                ///实例化KD树
                KDTree kD = new KDTree();
                ///最近K点
                List<pts> Points = new List<pts>();
                for (int i = 0; i < scatteredPositions.Length; i++)
                {
                    kPoints.Add(new pts((float)scatteredPositions[0, i], (float)scatteredPositions[1, i], (double)scatteredValues[0, i]));
                }
                for (int i = 0; i < k; i++)
                {
                    kD.CreateByPointList(kPoints);
                    pts middle = kD.FindNearest(new pts(target.X, target.Y, 0));
                    Points.Add(middle);
                    kPoints.Remove(middle);
                }

                ///在包含目标点之前重构插值点
                ///
                scatteredValues = new double[1, Points.Count];
                scatteredPositions = new double[2, Points.Count];
                for (int i = 0; i < Points.Count; i++)
                {
                    scatteredPositions[0, i] = Points[i].X;
                    scatteredPositions[1, i] = Points[i].Y;
                    scatteredValues[0, i] = Points[i].Value;
                }
                /*
                 * 包含目标点！
                 */
                Points.Add(new pts(target.X, target.Y, 0));
                /*
                 * 确定网格范围
                 */
                double start_x = Points[0].X;
                double end_x = Points[0].X;
                double start_y = Points[0].Y;
                double end_y = Points[0].Y;
                for (int i = 0; i < Points.Count; i++)
                {
                    start_x = Math.Min(Points[i].X, start_x);
                    end_x = Math.Max(Points[i].X, end_x);
                    start_y = Math.Min(Points[i].Y, start_y);
                    end_y = Math.Max(Points[i].Y, end_y);
                }
                start_x = Math.Round(start_x, retain);
                end_x = Math.Round(end_x, retain);
                start_y = Math.Round(start_y, retain);
                end_y = Math.Round(end_y, retain);
                ///构造函数
                Func<InArray<double>, InArray<double>, RetArray<double>> myFunc = (x, y) =>
                {
                    using (Scope.Enter(x, y))
                    {
                        return sin(x) * cos(y) * exp(-(x * x * y * y) / 4);
                    }
                };
                /*
                 * 创建网格
                 */
                Array<double> Y = 1, X = meshgrid(arange(start_x, step, end_x), arange(start_y, step, end_y), Y);
                Array<double> interpPositions = X[full].T.Concat(Y[full].T, 0);

                Array<double> Grid = reshape(X.C, X.S[0], X.S[1], 1);
                Grid[ellipsis, 1] = Y;
                Grid[ellipsis, 2] = myFunc(X, Y);
                /*
                 * 克里金插值模型
                 */
                Array<double> err = 1;
                Array<double> interpValues = Interpolation.kriging(scatteredValues, scatteredPositions, interpPositions, error: err);
                // let's plot it! 
                var scene = new Scene();
                var pc = scene.Add(new PlotCube(twoDMode: false));
                /*
                 * 绘制点
                 */
                Array<float> interpValues4Plot = tosingle(scatteredPositions.Concat(scatteredValues, 0));
                pc.Add(new ILNumerics.Drawing.Points()
                {
                    Positions = interpValues4Plot,
                    Color = Color.LightPink,
                    Size = 7
                });
                /*
                 * 将这些值恢复到合适的形状以进行绘图
                 */
                interpValues.a = reshape(interpValues, X.shape);
                err.a = reshape(err, X.shape);
                Array<double> outPut = X[":"].T.Concat(Y[":"].T, 0).Concat(interpValues[":"].T, 0);
                /*
                 * 绘制曲面
                 */
                pc.Add(new Surface(interpValues, X, Y, colormap: Colormaps.Jet, C: err)
                {
                    Fill = { Visible = false },
                    Wireframe = { Width = 2, Color = null }, // 'Color' controls the solid color. 'null' activates individual color mapped colors. 
                    Children = {
                    new Colorbar() {
                        // give the colorbar a title
                        new ILNumerics.Drawing.Label(text: @"Error, (\vartheta)") {
                            Anchor = new PointF(-.4f,1.1f)
                        }
                    }
                }
                });
                
                var cc = outPut.ToArray();
                for (long x_val = 0; x_val < cc.Length - 1; x_val++)
                {
                    var bb = Math.Round(cc[x_val], retain);
                    if ((Math.Round(cc[x_val], retain+1) == Math.Round(target.X, retain+1)) & (Math.Round(cc[x_val + 1], retain+1) == Math.Round(target.Y, retain+1)))
                    {
                        out_value = cc[x_val + 2];
                      //  MessageBox.Show("X:"+Math.Round(+cc[x_val], retain).ToString()+"Y:"+Math.Round(cc[x_val + 1], retain).ToString()+"vALUE:"+ out_value.ToString());
                    }
                }
                /*
                * 绘制目标点
                */
                Array<double> target_point = new double[2, 1] { { target.X }, { target.Y } };
                Array<double> target_value = new double[1, 1] { { out_value } };
                Array<float> interpValuest = tosingle(target_point.Concat(target_value, 0));
                pc.Add(new ILNumerics.Drawing.Points()
                {
                    Positions = interpValuest,
                    Color = Color.DarkRed,
                    Size = 9
                });
                // plot the original function (surface only, no wireframe)
                pc.Add(new Surface(Grid[ellipsis, 2], X, Y, colormap: Colormaps.Gray)
                {
                    Wireframe = { Visible = false, Color = Color.Black },
                    UseLighting = true
                });
                panel1.Scene = scene;
                // don't forget to prepare the scene for rendering + trigger a redraw!
                panel1.Configure();
                panel1.Refresh();
            }
            return out_value;
        }
        private void Form1_Load(object sender, EventArgs e)
        {

            Func<InArray<double>, InArray<double>, RetArray<double>> myFunc = (x, y) =>
            {
                return sin(x) * cos(y) * exp(-(x * x * y * y) / 4);
            };
            int ValuePoints = 50;
            Array<double> scatteredPositions = new double[2, ValuePoints];
            Array<double> scatteredValues = new double[1, ValuePoints];
            Random random = new Random(12);
            for (int i = 0; i < ValuePoints; i++)
            {
                scatteredPositions[0, i] = random.NextDouble() * 6 - 3;
                scatteredPositions[1, i] = random.NextDouble() * 6 - 3;
                scatteredValues[0, i] = myFunc(scatteredPositions[0, i], scatteredPositions[1, i]);
            }
          //  Kriging(scatteredValues, scatteredPositions, -3, 3, 0.1f);
           var cc = Kriging(scatteredValues, scatteredPositions, new PointF(0.132f, 1.2f), 5);
        }
    }

}
