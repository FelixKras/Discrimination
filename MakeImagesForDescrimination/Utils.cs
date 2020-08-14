using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using MakeImagesForDiscrimination;
using PixelFormat = System.Drawing.Imaging.PixelFormat;

namespace MakeImagesForDescrimination
{

    public class Rotations
    {
        /// <summary>
        /// Converts angle from degrees to radians and returns the result.
        /// </summary>
        /// <param name="angle">Angle in degrees to convert to radians.</param>
        /// <returns>Angle in radians.</returns>
        public static double ToRadians(double angle)
        {
            return angle * Math.PI / 180.0;
        }

        /// <summary>
        /// Returns 3D rotation matrix around Z axis with given angle 
        /// (expressed in radians).
        /// </summary>
        /// <param name="angle">Angle for rotation in radians.</param>
        /// <returns>3D rotation matrix around Z axis.</returns>

        public static double[,] CreateRotationMatrixZ(double angle)
        {
            return new double[3, 3]
            {
                {1.0, 0.0, 0.0},
                {0.0, Math.Cos(angle), -Math.Sin(angle)},
                {0.0, Math.Sin(angle), Math.Cos(angle)},

            };
        }


        public static double[,] CreateRotationPsiThetaPhi(double anglePsi, double angleThet, double anglePhi)
        {
            var cPsi = Math.Cos(anglePsi);
            var sPsi = Math.Sin(anglePsi);
            var cThet = Math.Cos(angleThet);
            var sThet = Math.Sin(angleThet);
            var cPhi = Math.Cos(anglePhi);
            var sPhi = Math.Sin(anglePhi);


            return new double[4, 4]
            {
                {cPsi*cThet,sPsi*cThet,-sThet,0},
                {cPsi*sThet*sPhi-sPsi*cPhi, sPsi*sThet*sPhi+cPsi*cPhi, cThet*sPhi,0},
                {cPsi*sThet*cPhi+sPsi*sPhi,sPsi*sThet*cPhi-cPsi*sPhi,cThet*cPhi,0},
                {0.0, 0.0,0.0,1}

            };
        }

        public static double[,] CreateRotationThetaPsiPhi(double angleThet, double anglePsi, double anglePhi)
        {
            var cPsi = Math.Cos(anglePsi);
            var sPsi = Math.Sin(anglePsi);
            var cThet = Math.Cos(angleThet);
            var sThet = Math.Sin(angleThet);
            var cPhi = Math.Cos(anglePhi);
            var sPhi = Math.Sin(anglePhi);

            return new double[4, 4]
            {
                {cPsi*cThet,-sPsi,cPsi*sThet,0},
                {-sPsi*cThet*cPhi+sThet*sPhi, -cPhi*cPsi,-sPsi*sThet*cPhi- cThet*sPhi,0},
                {sPsi*cThet*sPhi+sThet*cPhi,cPsi*sPhi,sPsi*sThet*sPhi-cThet*cPhi,0},
                {0.0, 0.0,0.0,1}

            };

            //return new double[4, 4]
            //{
            //    {cPsi*cThet,sPsi,-cPsi*sThet,0},
            //    {-sPsi*cThet*cPhi+sThet*sPhi, cPhi*cPsi,sPsi*sThet*cPhi+ cThet*sPhi,0},
            //    {sPsi*cThet*sPhi+sThet*cPhi,-cPsi*sPhi,-sPsi*sThet*sPhi+cThet*cPhi,0},
            //    {0.0, 0.0,0.0,1}
            //
            //};
        }

        //Rotation around X-axis
        public static double[,] CreateRotationMatrixRoll(double angle)
        {

            return new double[4, 4]
            {
                {1.0,0.0,0.0,0.0},
                {0.0,Math.Cos(angle), -Math.Sin(angle), 0.0},
                {0.0,Math.Sin(angle), Math.Cos(angle), 0.0},
                {0.0,0.0, 0.0, 1.0}
            };
        }

        //Rotation around Z-axis
        public static double[,] CreateRotationMatrixYaw(double angle)
        {
            return new double[4, 4]
            {
                {Math.Cos(angle), -Math.Sin(angle), 0.0, 0.0},
                {Math.Sin(angle), Math.Cos(angle), 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0},
                {0.0, 0.0, 0.0, 1.0}
            };
        }

        //Rotation around Y-axis
        public static double[,] CreateRotationMatrixPitch(double angle)
        {
            return new double[4, 4]
            {
                {Math.Cos(angle), 0.0, Math.Sin(angle),0.0},
                {0.0, 1.0, 0.0,0.0},
                {-Math.Sin(angle), 0.0, Math.Cos(angle),0.0},
                {0.0,0.0,0.0,1.0}
            };
        }

        //// <summary>
        //// Returns 3D rotation matrix around arbitrary unit vector defined by
        //// u with given angle (expressed in radians).
        //// </summary>
        //// <param name="angle">Angle for rotation in radians.</param>
        //// <returns>3D rotation matrix around general axis defined by u.</returns>
        //// <see cref="https://en.wikipedia.org/wiki/Rotation_matrix"/>
        public static double[,] CreateRotationMatrixGeneric(double[] u, double angle)
        {
            // Construct rotation matrix using formulas.
            // In simplified notation:
            // R = cos(angle) * I + sin(angle) * [u]_x + (1 - cos(angle)) * [u (x) u].
            // 
            // I is a 3x3 identity matrix.
            // 
            // [u]_x is a cross product matrix formed from u elements.
            // |0.0  -u_3  u_2|
            // |u_3  0.0  -u_1|
            // |-u_2 u_1   0.0|
            // 
            // u (x) u is a tensor product and equals uu^T.
            double a = Math.Cos(angle);
            double b = Math.Sin(angle);
            double c = 1 - a;

            double[,] part1 = new double[3, 3]
            {
                {a, 0.0, 0.0},
                {0.0, a, 0.0},
                {0.0, 0.0, a}
            };

            double[,] part2 = new double[3, 3]
            {
                {0.0, -b * u[2], b * u[1]},
                {b * u[2], 0.0, -b * u[0]},
                {-b * u[1], b * u[0], 0.0},
            };

            double[,] part3 = new double[3, 3]
            {
                {c * u[0] * u[0], c * u[0] * u[1], c * u[0] * u[2]},
                {c * u[1] * u[0], c * u[1] * u[1], c * u[1] * u[2]},
                {c * u[2] * u[0], c * u[2] * u[1], c * u[2] * u[2]}
            };

            double[,] R = BLAS.Add(BLAS.Add(part1, part2), part3);

            return R;
        }
    }

    public class Transltions
    {
        public static double[,] CreateTranslationMatrix(double[] trans_vector)
        {
            double[,] mat = new double[4, 4]
            {
                {1, 0.0, 0.0, trans_vector[0]},
                {0.0, 1, 0.0, trans_vector[1]},
                {0.0, 0.0, 1, trans_vector[2]},
                {0.0, 0.0, 0, 1.0}
            };
            return mat;
        }
    }
    public class Scales
    {
        public static double ScaleX;
        public static double ScaleY;
        public static double[,] ScaleMat;

        public static void SetScale(Size canvasSize)
        {
            ScaleX = (double)canvasSize.Width / Camera.CameraSettings.SensorWidth;
            ScaleY = (double)canvasSize.Height / Camera.CameraSettings.SensorHeight;
            ScaleMat = new double[4, 4]
            {
                {ScaleX, 0.0, 0.0, 0},
                {0.0, ScaleY, 0.0, 0},
                {0.0, 0.0, 1, 0},
                {0.0, 0.0, 0, 1.0}
            };

        }

    }
    public class Camera
    {
        public static double[,] CreateProjectionMatrix()
        {

            #region Explanation
            /*
            The intrinsic camera matrix is of the form:

            f_x s   x
            0   f_y y
            0   0   1
            
            Here, f_x and f_y are the focal lengths of the camera in the X and Y directions.
            s is the axis skew and is usually 0. x and y are the X and Y dimensions of the image
            produced by the camera, measured from the center of the image. 
            (So, they are half the length and width of the image.)
            
            We typically know the dimensions of the image produced by the camera.
            What is typically not provided are the focal lengths. Instead camera manufacturers provide
            the field of view (FOV) angle in the horizontal and vertical directions.

            Using the FOV angles, the focal lengths can be computed using trigonometry.
            For example, given the FOV a_x in the horizontal direction, 
            the focal length f_x can be computed using:

            f_x = x / tan(a_x / 2)

            We divide the FOV by 2 because this angle spans the entire horizontal or vertical view.

            As an example, consider the Primesense Carmine 1.09 depth camera. It produces a VGA (640×480) image. 
            Its specifications state a horizontal FOV of 57.5 degrees and vertical FOV of 45 degrees.

            Using the above information, we can compute its intrinsic camera matrix as:

            583.2829786373293       0.0                    320.0
            0.0                     579.4112549695428      240.0
            0.0                     0.0                     1.0
             
             */
            #endregion Explanation

            double[,] mat = new double[3, 4]
            {
                {CameraSettings.FocalX, 0.0, CameraSettings.SensorWidth/2D,0},
                {0.0, CameraSettings.FocalY, CameraSettings.SensorHeight/2D, 0},
                {0.0, 0.0, 1, 0}
                };
            return mat;
        }

        public static double[,] CreateAxisTranformMatrix()
        {
            return new double[,] { { 0, -1, 0, 0 }, { 0, 0, -1, 0 }, { 1, 0, 0, 0 }, { 0, 0, 0, 1 } };
        }
        public static class CameraSettings
        {
            // real fov angle is 32mRad =32D / 1000 * 180 / Math.PI; but we'll use much larger to simplify
            public static double Deg2Rad = Math.PI / 180;
            public static double FOVangEl = 0.1 / Deg2Rad;
            public static double FOVangAz = 0.1 / Deg2Rad;

            public static int SensorWidth = Logic.Settings.SensorSize.Width;
            public static int SensorHeight = Logic.Settings.SensorSize.Height;
            public static readonly double FocalX;
            public static readonly double FocalY;

            static CameraSettings()
            {
                FocalX = SensorWidth / 2D / Math.Tan(Rotations.ToRadians(FOVangAz) / 2);
                FocalY = SensorHeight / 2D / Math.Tan(Rotations.ToRadians(FOVangEl) / 2);
            }


        }

    }
    public class ConvexHull
    {
        [DllImport("kernel32.dll")]
        private static extern IntPtr LoadLibrary(string sDllToLoad);
        [DllImport("kernel32.dll")]
        private static extern bool FreeLibrary(IntPtr hModule);

        [DllImport("concaveman.dll", EntryPoint = "pyconcaveman2d", CallingConvention = CallingConvention.Cdecl)]
        public static extern void pyconcaveman2d([In, Out] strPoint[] points_c, uint num_points,
            [In, Out] int[] hull_points_c, [In, Out] uint num_hull_points,
            double concavity, double lengthThreshold,
            [In, Out]  strPoint[] concave_points_c, ref uint num_concave_points);


        private IntPtr hUnMannagedLib;
        public void Dispose()
        {
            FreeLibrary(hUnMannagedLib);
        }

        ~ConvexHull()
        {
            FreeLibrary(hUnMannagedLib);
        }

        private static IntPtr LoadApi()
        {
            IntPtr handle;
            try
            {
                string sDllName = @"concaveman.dll";
                handle = LoadLibrary(sDllName);
            }
            catch (Exception ex)
            {
                handle = IntPtr.Zero;
                throw;
            }
            return handle;
        }

        public static List<PointOnGrid> GetConcaveHull(List<PointOnGrid> points)
        {
            /*
            strPoint[] inPoints = points.Cnvrt2strPnts();

            List<IndexedPoint> convexHullPnts = ConvexHull.GetConvexHull(inPoints.Cnvrt2LstIndxPnts());

            List<PointF> convexHullPntsF = ConvexHull.GetConvexHull(pointf);

            int[] hullPoints = convexHullPnts.Select(point => point.Index).ToArray();

            strPoint[] outPoints = new strPoint[inPoints.Length];

            uint numHullPts = (uint)hullPoints.Length;
            uint numInputPoints = (uint)inPoints.Length;
            uint numOutPts = (uint)inPoints.Length;
            double concavity = 1;
            double lengthThresh = 1;//(uint)inPoints.Length;
            

            pyconcaveman2d(inPoints, numInputPoints, hullPoints, numHullPts,
                concavity, lengthThresh, outPoints, ref numOutPts);

            return outPoints;
            */
            return null;
        }

        private static double cross(PointF O, PointF A, PointF B)
        {
            return (A.X - O.X) * (B.Y - O.Y) - (A.Y - O.Y) * (B.X - O.X);
        }
        public static List<PointF> GetConvexHull(List<PointF> points)
        {
            if (points == null)
                return null;

            if (points.Count() <= 1)
                return points;

            int n = points.Count(), k = 0;
            List<PointF> H = new List<PointF>(new PointF[2 * n]);

            points.Sort((a, b) =>
            {
                if (Math.Abs(a.X - b.X) < 1e-6)
                    return a.Y.CompareTo(b.Y);
                else
                    return a.X.CompareTo(b.X);
            });

            // Build lower hull
            for (int i = 0; i < n; ++i)
            {
                while (k >= 2 && cross(H[k - 2], H[k - 1], points[i]) <= 0)
                    k--;
                H[k++] = points[i];
            }

            // Build upper hull
            for (int i = n - 2, t = k + 1; i >= 0; i--)
            {
                while (k >= t && cross(H[k - 2], H[k - 1], points[i]) <= 0)
                    k--;
                H[k++] = points[i];
            }

            return H.Take(k - 1).ToList();
        }

        private static double cross(PointOnGrid O, PointOnGrid A, PointOnGrid B)
        {
            return (A.X - O.X) * (B.Y - O.Y) - (A.Y - O.Y) * (B.X - O.X);
        }
        public static List<PointOnGrid> GetConvexHull(List<PointOnGrid> points)
        {
            if (points == null)
                return null;

            if (points.Count() <= 1)
                return points;

            int n = points.Count(), k = 0;
            List<PointOnGrid> H = new List<PointOnGrid>(new PointOnGrid[2 * n]);

            points.Sort((a, b) =>
                a.X == b.X ? a.Y.CompareTo(b.Y) : a.X.CompareTo(b.X));

            // Build lower hull
            for (int i = 0; i < n; ++i)
            {
                while (k >= 2 && cross(H[k - 2], H[k - 1], points[i]) <= 0)
                    k--;
                H[k++] = points[i];
            }

            // Build upper hull
            for (int i = n - 2, t = k + 1; i >= 0; i--)
            {
                while (k >= t && cross(H[k - 2], H[k - 1], points[i]) <= 0)
                    k--;
                H[k++] = points[i];
            }

            return H.Take(k - 1).ToList();
        }

        public static List<int> GetLocalConvex(List<PointOnGrid> points, double scale)
        {

            #region Navot Israeli code
            /*
            %computes the envelope of a set of points in the x,y plane
            %uses a variation of convex hull ahere at each step the next leg is taken
            %as the "local" continuation of the convex hull. Namely, the continuation
            %of convhull with the points in the local neighborhood of the current node
            %
            % input:
            % x,y are points coordinates
            % scale is the radius of the local neighborhood
            % output:
            % env_ii is an array of indices into x,y. Namely, x(env_ii), y(env_ii) is the
            % envelope polygon
            
            function env_ii=local_convhull_envelope(x,y,scale)
            
            
            [~, start_p]=min(x);
            
            p=start_p;
            in_dir=[0; 1];
            
            env_x=x(p);
            env_y=y(p);
            env_ii=p;
            
            
            while (1)
                d=sqrt(((x-x(p))/scale).^2+((y-y(p))/scale).^2);
                ii=find(d<=1);
                
                
                min_alpha=1000;
                best_q=-1;
                best_out_dir=[];
                in_dir_n=[in_dir(2); -in_dir(1)];
                
                for qii=1:length(ii)
                    q=ii(qii);
                    if q~=p
                        out_dir=[x(q)-x(p); y(q)-y(p)];
                        out_dir=out_dir/norm(out_dir);
                        
                        alpha=atan2(dot(out_dir,in_dir_n),dot(out_dir,in_dir));
                        if (alpha<min_alpha)
                            intersects=0;
                            for k=1:length(env_x)-2
                                A=[(env_x(k+1)-env_x(k)) (x(p)-x(q)); (env_y(k+1)-env_y(k)) (y(p)-y(q))];
                                if (det(A)>1e-16)
                                    V=[x(p)-env_x(k); y(p)-env_y(k)];
                                    s=A\V;
                                    
                                    if ((s(1)>0) && (s(1)<1) && (s(2)>0) && (s(2)<1))
                                        intersects=1;
                                        break;
                                    end
                                end
                            end
                            if intersects==0
                                min_alpha=alpha;
                                best_q=q;
                                best_out_dir=out_dir;
                            end
                        end
                    end
                end
                
                env_x(end+1)=x(best_q);
                env_y(end+1)=y(best_q);
                env_ii(end+1)=best_q;
                p=best_q;
                in_dir=best_out_dir;
                
                
                
                if p==start_p;
                    break;
                end
            end
            
            
            
            
            
             */
            #endregion Navot Israeli code

            int startP = points.IndexOfMin();                        //[~, start_p]=min(x);
            int pInd = startP;                                                         //p=start_p;
            double[] inDir = new double[] { 0, 1 };                                  //in_dir=[0; 1];

            List<double> env_x = new List<double>() { points[pInd].X };                         //env_x=x(p);
            List<double> env_y = new List<double>() { points[pInd].Y };                         //env_y=y(p);
            List<int> env_ii = new List<int>() { pInd };                                  //env_ii=p;

            int iter = 0; //debug
            int best_q;
            double min_alpha;
            double alpha;

            List<int> closeIndexes = new List<int>();
            List<double> dist;

            double[] best_out_dir = new double[2];
            double[] in_dir_n = new double[2] { inDir[1], -inDir[0] };
            double[] outDir = new double[2];
            double[] V = new double[2];
            double[] s = new double[2];
            double[,] A = new double[2, 2];
            while (true)
            {
                iter++;
                if (iter > 1000 && iter > points.Count * 0.3)
                {
                    //env_ii.Clear();
                    break;
                }
                dist = points.Select(p => Math.Sqrt(Math.Pow((p.X - points[pInd].X) / scale, 2) +
                                                                      Math.Pow((p.Y - points[pInd].Y) / scale, 2))).ToList();    //d=sqrt(((x-x(p))/scale).^2+((y-y(p))/scale).^2);


                dist.FindIndexes(f => f <= 1, ref closeIndexes);//ii = find(d <= 1);


                min_alpha = 1000;                                   //min_alpha=1000;
                best_q = -1;                                        //best_q=-1;

                //best_out_dir = new double[2];                           //best_out_dir=[];
                best_out_dir[0] = 0;
                best_out_dir[1] = 0;
                in_dir_n[0] = inDir[1];
                in_dir_n[1] = -inDir[0];        //in_dir_n=[in_dir(2); -in_dir(1)];
                outDir[0] = 0;
                outDir[1] = 0;
                for (int qi = 0; qi < closeIndexes.Count; qi++)
                {
                    int q = closeIndexes[qi]; //q=ii(qii);
                    if (q != pInd)
                    {
                        //out_dir=[x(q)-x(p); y(q)-y(p)];
                        outDir[0] = points[q].X - points[pInd].X;
                        outDir[1] = points[q].Y - points[pInd].Y;

                        outDir.DivideByNorm(); //out_dir=out_dir/norm(out_dir);
                        alpha = Math.Atan2(dot(outDir, in_dir_n), dot(outDir, inDir)); //alpha=atan2(dot(out_dir,in_dir_n),dot(out_dir,in_dir));

                        if (alpha < min_alpha) //if (alpha<min_alpha)
                        {
                            int intersects = 0;  //    intersects=0;
                            for (int kk = 0; kk <= env_x.Count - 3; kk++)//    for k=1:length(env_x)-2
                            {
                                //double[,] A = new double[,] { { env_x[kk + 1] - env_x[kk], points[pInd].X - points[q].X }, { env_y[kk + 1] - env_y[kk], points[pInd].Y - points[q].Y } }; //        A=[(env_x(k+1)-env_x(k)) (x(p)-x(q)); (env_y(k+1)-env_y(k)) (y(p)-y(q))];
                                A[0, 0] = env_x[kk + 1] - env_x[kk]; // A=[(env_x(k+1)-env_x(k)) (x(p)-x(q)); (env_y(k+1)-env_y(k)) (y(p)-y(q))];
                                A[0, 1] = points[pInd].X - points[q].X;
                                A[1, 0] = env_y[kk + 1] - env_y[kk];
                                A[1, 1] = points[pInd].Y - points[q].Y;

                                if (det(A) > 1e-16) //if (det(A)>1e-16)
                                {
                                    V[0] = points[pInd].X - env_x[kk];  //V =[x(p) - env_x(k); y(p) - env_y(k)];
                                    V[1] = points[pInd].Y - env_y[kk];  //V =[x(p) - env_x(k); y(p) - env_y(k)];
                                    Pinv(A, V, ref s);  //s=A\V;

                                    if ((s[0] > 0) && (s[0] < 1) && (s[1] > 0) && (s[1] < 1))      //if ((s(1)>0) && (s(1)<1) && (s(2)>0) && (s(2)<1))
                                    {
                                        intersects = 1;        // intersects=1;
                                        break;                 // break;
                                    }
                                }
                            }

                            if (intersects == 0)
                            {
                                min_alpha = alpha; // min_alpha = alpha;
                                best_q = q;        // best_q=q;
                                best_out_dir[0] = outDir[0]; //        best_out_dir=out_dir;
                                best_out_dir[1] = outDir[1];
                            }
                        }

                    }

                }
                env_x.Add(points[best_q].X); //env_x(end + 1) = x(best_q);
                env_y.Add(points[best_q].Y); //env_y(end + 1) = y(best_q);
                env_ii.Add(best_q);          //env_ii(end + 1) = best_q;
                pInd = best_q;           //p = best_q;
                inDir[0] = best_out_dir[0];        //in_dir = best_out_dir;
                inDir[1] = best_out_dir[1];

                if (pInd == startP) //if p == start_p;
                {
                    break;
                }
            }//end
            return env_ii;
        }

        private static void Pinv(double[,] A, double[] b, ref double[] s)
        {
            s[0] = (A[1, 1] * b[0] - A[0, 1] * b[1]) / (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]);
            s[1] = (A[0, 0] * b[1] - A[1, 0] * b[0]) / (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]);
        }

        public static double dot(double[] a, double[] b)
        {
            double result = 0;
            for (int ii = 0; ii < a.Length; ii++)
            {
                result += a[ii] * b[ii];
            }
            return result;
        }
        private static float Norm(float[] inVec)
        {
            double result = 0;
            for (int ii = 0; ii < inVec.Length; ii++)
            {
                result += inVec[ii] * inVec[ii];
            }
            return (float)Math.Sqrt(result);
        }

        public static double det(double[,] A)
        {
            double result = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0];
            return result;
        }
    }

    public class BLAS
    {
        /// <summary>
        /// Adds two vectors of same length and returns result (a + b).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>a + b.</returns>
        public static double[] Add(double[] a, double[] b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        /// <summary>
        /// Adds two vectors of same length and returns result (a + b).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>a + b.</returns>
        public static int[] Add(int[] a, int[] b)
        {
            int[] result = new int[a.Length];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        /// <summary>
        /// Adds two matrices of same size and returns result (A + B).
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns>A + B.</returns>
        public static double[,] Add(double[,] A, double[,] B)
        {
            double[,] result = new double[A.GetLength(0), A.GetLength(1)];

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    result[i, j] = A[i, j] + B[i, j];
                }
            }

            return result;
        }

        /// <summary>
        /// Subtracts two vectors of same length and returns result (a - b).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>a - b.</returns>
        public static double[] Subtract(double[] a, double[] b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = a[i] - b[i];
            }

            return result;
        }

        /// <summary>
        /// Subtracts two vectors of same length and returns result (a - b).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>a - b.</returns>
        public static int[] Subtract(int[] a, int[] b)
        {
            int[] result = new int[a.Length];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = a[i] - b[i];
            }

            return result;
        }

        /// <summary>
        /// Returns the L2 norm of a vector ||v||.
        /// </summary>
        /// <param name="v">Vector to compute L2-norm of.</param>
        /// <returns>L2-norm of vector ||v||.</returns>
        public static double Norm2(double[] v)
        {
            double result2 = 0.0;

            foreach (double item in v)
            {
                result2 += item * item;
            }

            return Math.Sqrt(result2);
        }

        /// <summary>
        /// Returns the L2 norm of a vector ||v||.
        /// </summary>
        /// <param name="v">Vector to compute L2-norm of.</param>
        /// <returns>L2-norm of vector ||v||.</returns>
        public static int Norm2(int[] v)
        {
            int result2 = 0;

            foreach (int item in v)
            {
                result2 += item * item;
            }

            return (int)Math.Sqrt(result2);
        }

        /// <summary>
        /// In-place normalization of given vector to become unit-length.
        /// </summary>
        /// <param name="v">Vector to normalize.</param>
        public static void Normalize(double[] v)
        {
            double norm = Norm2(v);

            for (int i = 0; i < v.Length; i++)
            {
                v[i] /= norm;
            }
        }

        /// <summary>
        /// Returns the dot product between two vectors (a, b).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>Dot product (a, b).</returns>
        public static double Dot(double[] a, double[] b)
        {
            double result = 0.0;

            for (int i = 0; i < a.Length; i++)
            {
                result += a[i] * b[i];
            }

            return result;
        }

        /// <summary>
        /// GEMM operation, matrix-matrix multiplication (A * B).
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns>A * B.</returns>

        public static double[,] Multiply(double[,] A, double[,] B)
        {
            double[,] result = new double[A.GetLength(0), B.GetLength(1)];

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < B.GetLength(1); j++)
                {
                    for (int k = 0; k < A.GetLength(1); k++)
                    {
                        result[i, j] += A[i, k] * B[k, j];
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// GEMV operation, matrix-vector multiplication (A * b).
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns>A * b.</returns>
        public static double[] Multiply(double[,] A, double[] b)
        {
            double[] result = new double[b.Length];

            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < b.Length; j++)
                {
                    result[i] += A[i, j] * b[j];
                }
            }

            return result;
        }

        /// <summary>
        /// GEMV operation, vector-matrix multiplication (a * B).
        /// </summary>
        /// <param name="a"></param>
        /// <param name="B"></param>
        /// <returns>a * B.</returns>
        public static double[] Multiply(double[] a, double[,] B)
        {
            double[] result = new double[B.GetLength(0)];

            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < a.Length; j++)
                {
                    result[i] += a[j] * B[j, i];
                }
            }

            return result;
        }

        /// <summary>
        /// Multiplies each vector element with scalar and returns the result (a .* b).
        /// </summary>
        /// <param name="a">Vector.</param>
        /// <param name="b">Scalar.</param>
        /// <returns>a .* b.</returns>
        public static double[] PointwiseMultiply(double[] a, double b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] * b;
            }

            return result;
        }

        /// <summary>
        /// Divides each vector element with scalar and returns the result (a ./ b).
        /// </summary>
        /// <param name="a">Vector.</param>
        /// <param name="b">Scalar.</param>
        /// <returns>a ./ b.</returns>
        public static double[] PointwiseDivide(double[] a, double b)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] / b;
            }

            return result;
        }

        /// <summary>
        /// Computes the absolute value of each vector element and returns the result.
        /// </summary>
        /// <param name="a"></param>
        /// <returns>Absolute value of each vector element.</returns>
        public static double[] PointwiseAbsolute(double[] a)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = Math.Abs(a[i]);
            }

            return result;
        }

        /// <summary>
        /// Rounds each vector element and returns the result.
        /// </summary>
        /// <param name="a"></param>
        /// <returns>Rounded value of each vector element.</returns>
        public static double[] PointwiseRound(double[] a)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = Math.Round(a[i]);
            }

            return result;
        }

        /// <summary>
        /// Compute floor for each vector element and returns the result.
        /// </summary>
        /// <param name="a"></param>
        /// <returns>Floored value of each vector element.</returns>
        public static double[] PointwiseFloor(double[] a)
        {
            double[] result = new double[a.Length];

            for (int i = 0; i < a.Length; i++)
            {
                result[i] = Math.Floor(a[i]);
            }

            return result;
        }
    }

    public static class Utils
    {
        public static int IndexOfMin(this List<PointOnGrid> sequence)
        {
            int index = 0;
            float MinValue = float.MaxValue;
            for (int ii = 0; ii < sequence.Count(); ii++)
            {
                if (sequence[ii].X < MinValue)
                {
                    MinValue = sequence[ii].X;
                    index = ii;
                }
            }

            return index;
        }
        public static BitmapSource ArrayToGrayBitmap(ushort[,] SourceImage)
        {
            const int icColorCount = 2;
            BitmapSource source;
            int arrayWidth = SourceImage.GetLength(1), arrayHeight = SourceImage.GetLength(0);

            int iPaddedWidth = arrayWidth;
            if (arrayWidth % 4 != 0)
            {
                iPaddedWidth = (arrayWidth / 4 + 1) * 4;
            }

            Bitmap bmp = new Bitmap(arrayWidth, arrayHeight, PixelFormat.Format16bppGrayScale);
            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadWrite, bmp.PixelFormat);
            int size = bitmapData.Stride * bitmapData.Height;
            byte[] DataForBitmapOutBytes = new byte[size];
            int col, row, kk;
            try
            {
                // stride is padded to multiple of 4. iterate over the source not destination, otherwise out of range exception

                for (int ii = 0; ii < SourceImage.Length; ii++)
                {
                    col = ii % arrayWidth;
                    row = ii / arrayWidth;
                    if (row > 0)
                    {

                    }
                    kk = row * iPaddedWidth * icColorCount + col * icColorCount;
                    DataForBitmapOutBytes[kk] = (byte)(SourceImage[row, col]);
                    DataForBitmapOutBytes[kk + 1] = (byte)((SourceImage[row, col] & 0xff00) >> 8);


                }
                Marshal.Copy(DataForBitmapOutBytes, 0, bitmapData.Scan0, DataForBitmapOutBytes.Length);
            }
            catch (Exception e)
            {

                throw;
            }
            finally
            {
                /* pixel formats
                 switch (bmpPixelFormat)
                 {
                 case PixelFormat.Format32bppArgb:
                 format = PixelFormats.Bgr32;
                 break;
                 case PixelFormat.Format8bppIndexed:
                 format = PixelFormats.Gray8;
                 break;
                 case PixelFormat.Format16bppGrayScale:
                 format = PixelFormats.Gray16;
                 break;
                 }
                */
                source = BitmapSource.Create(bmp.Width, bmp.Height,
                   bmp.HorizontalResolution,
                   bmp.VerticalResolution,
                   PixelFormats.Gray16
                   , null,
                   bitmapData.Scan0,
                   bitmapData.Stride * bmp.Height,
                   bitmapData.Stride);

                bmp.UnlockBits(bitmapData);
            }


            return source;
        }

        public static Bitmap ArrayToBitmap(ushort[,] SourceImage, int[][] palette)
        {
            const int icColorCount = 3;
            int arrayWidth = SourceImage.GetLength(1), arrayHeight = SourceImage.GetLength(0);

            int iPaddedWidth = arrayWidth;
            if (arrayWidth % 4 != 0)
            {
                iPaddedWidth = (arrayWidth / 4 + 1) * 4;
            }

            Bitmap bmp = new Bitmap(arrayWidth, arrayHeight, PixelFormat.Format24bppRgb);
            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadWrite, bmp.PixelFormat);
            int bytesPerPixel = Image.GetPixelFormatSize(bitmapData.PixelFormat) / 8;
            int size = bitmapData.Stride * bitmapData.Height;
            byte[] DataForBitmapOutBytes = new byte[size];
            int col, row;
            try
            {
                // stride is padded to multiple of 4. iterate over the source not destination, otherwise out of range exception

                for (int ii = 0; ii < SourceImage.Length; ii++)
                {
                    col = ii % arrayWidth;
                    row = ii / arrayWidth;
                    if (row > 0)
                    {

                    }
                    int kk = row * iPaddedWidth * icColorCount + col * icColorCount;
                    DataForBitmapOutBytes[kk] = (byte)palette[SourceImage[row, col]][0];
                    DataForBitmapOutBytes[kk + 1] = (byte)palette[SourceImage[row, col]][1];
                    DataForBitmapOutBytes[kk + 2] = (byte)palette[SourceImage[row, col]][2];

                }
                Marshal.Copy(DataForBitmapOutBytes, 0, bitmapData.Scan0, DataForBitmapOutBytes.Length);
            }
            catch (Exception e)
            {

                throw;
            }
            finally
            {
                bmp.UnlockBits(bitmapData);
            }
            return bmp;
        }
        public static int[][] LoadPalette(string i_FileName)
        {
            const int icColorCount = 3;
            string[] paletteStringArray = File.ReadAllLines(i_FileName, Encoding.ASCII);
            int[][] paletteResult = new int[paletteStringArray.Length][];
            int R, G, B;

            try
            {
                for (int ii = 0; ii < paletteStringArray.Length; ii++)
                {
                    string[] RGB = paletteStringArray[ii].Split(new char[] { ',', ' ' });
                    R = int.Parse(RGB[1]);
                    G = int.Parse(RGB[2]);
                    B = int.Parse(RGB[3]);
                    paletteResult[ii] = new int[] { B, G, R };
                }
            }
            catch (Exception e)
            {
                throw;
            }
            return paletteResult;
        }

        public static T RandomIndex<T>(this T[] arr)
        {
            Random rnd = new Random(DateTime.Now.Millisecond);
            int length = arr.Length;
            int index = rnd.Next(0, length);
            return arr[index];

        }

        public static void SaveTo16bPNG(this BitmapSource bmpSourceGray, string spath)
        {
            using (FileStream fs = new FileStream(spath, FileMode.Create, FileAccess.Write))
            {
                PngBitmapEncoder encoder = new PngBitmapEncoder();
                //BitmapMetadata bmpMeta=new BitmapMetadata("png");
                //bmpMeta.SetQuery("/Text/Description", description);
                encoder.Frames.Add(BitmapFrame.Create(bmpSourceGray));
                encoder.Save(fs);
            }
        }

        public static void WriteCSV(this string[] content, string fname)
        {
            using (FileStream fs = new FileStream(fname, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite))
            using (StreamWriter sw = new StreamWriter(fs))
            {
                foreach (var line in content)
                {
                    sw.WriteLine(line);
                }
            }
        }
        public static string[] ToCSV(this List<PointOnGrid> lst)
        {

            string[] result = new string[lst.Count + 1];
            result[0] = "Index, X,Y,Component";
            for (int ii = 0; ii < lst.Count; ii++)
            {
                result[ii + 1] = lst[ii].Index + "," + lst[ii].X + "," + lst[ii].Y + "," + lst[ii].Component;
            }

            return result;
        }

        public static List<PointF> Cnvrt2LstPntf(this strPoint[] stPoints)
        {
            List<PointF> result = new List<PointF>(stPoints.Length);
            for (int ii = 0; ii < stPoints.Length; ii++)
            {
                result.Add(new PointF((float)stPoints[ii].x, (float)stPoints[ii].y));
            }

            return result;
        }

        public static List<PointOnGrid> Cnvrt2LstIndxPts(this strPoint[] stPoints)
        {
            List<PointOnGrid> result = new List<PointOnGrid>(stPoints.Length);
            for (int ii = 0; ii < stPoints.Length; ii++)
            {
                result.Add(new PointOnGrid((float)stPoints[ii].x, (float)stPoints[ii].y, 0, ii));
            }

            return result;
        }

        public static strPoint[] Cnvrt2strPnts(this List<PointOnGrid> Points)
        {
            strPoint[] result = new strPoint[Points.Count];
            for (int ii = 0; ii < Points.Count; ii++)
            {
                result[ii].x = Points[ii].X;
                result[ii].y = Points[ii].Y;
            }

            return result;
        }

        public static void FindIndexes<T>(this List<T> Elements, Func<T, bool> func, ref List<T> Indexes)
        {
            Indexes.Clear();
            for (int ii = 0; ii < Indexes.Count; ii++)
            {
                if (func(Elements[ii]))
                {
                    Indexes.Add(Elements[ii]);
                }
            }
        }
        public static void FindIndexes<T>(this List<T> Floats, Func<T, bool> func, ref List<int> closeIndexes)
        {
            closeIndexes.Clear();
            for (int ii = 0; ii < Floats.Count; ii++)
            {
                if (func(Floats[ii]))
                {
                    closeIndexes.Add(ii);
                }
            }
        }

        public static void DivideByNorm(this double[] a)
        {
            double result = 0;
            for (int ii = 0; ii < a.Length; ii++)
            {
                result += a[ii] * a[ii];
            }
            result = Math.Sqrt(result);
            for (int ii = 0; ii < a.Length; ii++)
            {
                a[ii] /= result;
            }
        }

        public static void CreateImage(this List<PointOnGrid> points, string spath)
        {
            const int icColorCount = 3;


            int maxXShape = (int)points.Max(p => p.X);
            int maxYShape = (int)points.Max(p => p.Y);
            int minXShape = (int)points.Min(p => p.X);
            int minYShape = (int)points.Min(p => p.Y);

            byte[,] SourceImage = new byte[maxXShape - minXShape + 1, maxYShape - minYShape + 1];

            for (int ii = 0; ii < points.Count; ii++)
            {
                SourceImage[((int)points[ii].X - minXShape), ((int)points[ii].Y - minYShape)] = 250;
            }
            int arrayWidth = SourceImage.GetLength(1);
            int arrayHeight = SourceImage.GetLength(0);

            int iPaddedWidth = arrayWidth;
            if (arrayWidth % 4 != 0)
            {
                iPaddedWidth = (arrayWidth / 4 + 1) * 4;
            }

            Bitmap bmp = new Bitmap(arrayWidth, arrayHeight, PixelFormat.Format24bppRgb);
            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadWrite, bmp.PixelFormat);
            int bytesPerPixel = Image.GetPixelFormatSize(bitmapData.PixelFormat) / 8;
            int size = iPaddedWidth * icColorCount * bitmapData.Height;
            byte[] DataForBitmapOutBytes = new byte[size];
            int col, row;
            try
            {
                // stride is padded to multiple of 4. iterate over the source not destination, otherwise out of range exception

                for (int ii = 0; ii < SourceImage.Length - 3; ii++)
                {
                    col = ii % arrayWidth;
                    row = ii / arrayWidth;


                    int kk = row * iPaddedWidth * icColorCount + col * icColorCount;

                    DataForBitmapOutBytes[kk] = SourceImage[row, col];
                    DataForBitmapOutBytes[kk + 1] = SourceImage[row, col];
                    DataForBitmapOutBytes[kk + 2] = SourceImage[row, col];
                }
                Marshal.Copy(DataForBitmapOutBytes, 0, bitmapData.Scan0, DataForBitmapOutBytes.Length);
            }
            catch (Exception e)
            {

                throw;
            }
            finally
            {
                bmp.UnlockBits(bitmapData);
            }
            bmp.Save(spath);

        }


    }
}
