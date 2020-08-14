using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using MakeImagesForDescrimination;
using PixelFormat = System.Drawing.Imaging.PixelFormat;

namespace MakeImagesForDiscrimination
{
    class Cylinder : baseShape, IDrawable
    {
        public Cylinder(float H, float W)
        {
            Vertices = new List<PointOnGrid>();
            Vertices.Add(new PointOnGrid(-W / 2, -H / 2, 0, 5));
            Vertices.Add(new PointOnGrid(W / 2, -H / 2, 0, 5));
            Vertices.Add(new PointOnGrid(W / 2, H / 2, 0, 10));
            Vertices.Add(new PointOnGrid(-W / 2, H / 2, 0, 10));

            oMaxMin = Logic.Settings.Temperatures.EngTemp;

        }

        public override uint GetShapeID()
        {
            return 2;
        }
    }
    class Conus : baseShape, IDrawable
    {
        public Conus(float H, float W)
        {
            Vertices = new List<PointOnGrid>();
            Vertices.Add(new PointOnGrid(0, H * 2 / 3, 0, 10));
            Vertices.Add(new PointOnGrid(-W / 2, -H * 1 / 3, 0, 5));
            Vertices.Add(new PointOnGrid(W / 2, -H * 1 / 3, 0, 5));

            oMaxMin = Logic.Settings.Temperatures.RvTemp;
        }

        public override uint GetShapeID()
        {
            return 1;
        }
    }
    class ComplexShape : baseShape, IDrawable
    {
        public static bool IsWriteDebug = false;
        private class StateVector
        {
            internal double RollAngle;
            internal double YawAngle;
            internal double PitchAngle;
            internal double XCG;
            internal double YCG;
            internal double ZCG;
        }
        private uint ShapeID;
        StateVector stateVec;
        public ComplexShape(List<PointOnGrid> inputPoints, MaxMin maxmin, eShape shape)
        {
            Vertices = inputPoints;
            //oMaxMin = new MaxMin(10000, 8000);
            oMaxMin = maxmin;
            stateVec = new StateVector()
            { PitchAngle = 0, RollAngle = 0, YawAngle = 0, XCG = 100, YCG = 0, ZCG = 0 };
            ShapeID = (uint)shape;

        }

        public override void Render(ref ushort[,] canvas, ref RectangleF rect)
        {
            //double area = 0;
            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);


            List<PointOnGrid> projectedPoints = ComputeProjection(stateVec, Vertices);
            if (IsWriteDebug)
                File.WriteAllLines(
                    @"C:\Users\88felix\source\repos\MakeImagesForDescrimination\TempDebug\projectedPoints.csv",
                    projectedPoints.ToCSV());
            
            if (projectedPoints.Count < 0.7 * Vertices.Count)
            {
                rect = RectangleF.Empty;
            }
            else
            {
                maxXShape = projectedPoints.Max(p => p.X);
                maxYShape = projectedPoints.Max(p => p.Y);
                minXShape = projectedPoints.Min(p => p.X);
                minYShape = projectedPoints.Min(p => p.Y);

                if (Logic.Settings.UseConcHull)
                {
                    double len = Math.Sqrt((maxXShape - minXShape) * (maxXShape - minXShape) +
                                        (maxYShape - minYShape) * (maxYShape - minYShape));
                    List<int> newCHull = ConvexHull.GetLocalConvex(projectedPoints, len / 10);
                    ushort[,] tempArray = new ushort[227, 227];

                    #region Obsolete, get convex hull by parts

                    // get convex hull by parts
                    //List<List<PointOnGrid>> convexHulls = new List<List<PointOnGrid>>();
                    //for (int ii = 1; ii < 10; ii++)
                    //{
                    //List<PointOnGrid> lstPart = (projectedPoints.Where(p => p.Component == ii)).ToList();
                    //if (lstPart.Count > 0)
                    //{
                    //    List<PointOnGrid> partialContour = ConvexHull.GetConvexHull(lstPart);
                    //    if (IsWriteDebug)
                    //        File.WriteAllLines(@"C:\Users\88felix\source\repos\MakeImagesForDescrimination\TempDebug\partialConv" + ii + ".csv", partialContour.ToCSV());

                    //    //PlotPoints(partialContour,ref tempArray);

                    //    convexHulls.Add(partialContour);
                    //}
                    //}
                    //List<PointOnGrid> mergedHull = MergeConvexHulls(convexHulls);

                    #endregion Obsolete, get convex hull by parts

                    List<PointOnGrid> mergedHull = GetPointsByIndex(projectedPoints, newCHull);
                    //mergedHull.CreateImage("c:\\temp\\temp.bmp");

                    //ShiftFromEdge(Vertices, canvas);


                    ScaleTargetAndAddToCanvas(canvas, mergedHull, oMaxMin, ref rect);

                }

                #region Obsolete, use cpp version of concave hull

                if (false)
                {
                    //List<PointOnGrid> reducedPoints = projectedPointsContour.Where((x, i) => i % 20 == 0).ToList();
                    //int[] hullPointsIndxs = reducedPoints.Select(point => point.Index).ToArray();
                    //strPoint[] outPoints = new strPoint[reducedPoints.Count];
                    //uint numHullPts = (uint)hullPointsIndxs.Length;
                    //uint numInputPoints = (uint)reducedPoints.Count;
                    //uint numOutPts = (uint)reducedPoints.Count;
                    //double concavity = 1;
                    //double lengthThresh = 1;//(uint)inPoints.Length;
                    //ConvexHull.pyconcaveman2d(projectedPoints.Cnvrt2strPnts(), numInputPoints, hullPointsIndxs, numHullPts,
                    //    concavity, lengthThresh, outPoints, ref numOutPts);

                }

                #endregion Obsolete, get convex hull by parts

                //}
                if(false) //fill by dense points
                {
                    ScaleTargetAndAddToCanvasNew(canvas, projectedPoints, oMaxMin, ref rect);
                }




            }
        }

        internal void ScaleTargetAndAddToCanvasNew(ushort[,] canvas, List<PointOnGrid> projectedPoints, MaxMin maxMin, ref RectangleF rect)
        {
            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);

            double maxX = ValidateSize(projectedPoints.Max(p => p.X), width);
            double maxY = ValidateSize(projectedPoints.Max(p => p.Y), height);
            double minX = ValidateSize(projectedPoints.Min(p => p.X), width);
            double minY = ValidateSize(projectedPoints.Min(p => p.Y), height);

            rect.X = (float)minX;
            rect.Y = (float)minY;
            rect.Width = (float)(maxX - minX);
            rect.Height = (float)(maxY - minY);


            var xcg = minX + (maxX - minX) / 2;
            var ycg = minY + (maxY - minY) / 2;
            var RadiusAroundShape = Math.Sqrt(((maxX - minX) / 2) * ((maxX - minX) / 2) +
                                             ((maxY - minY) / 2) * ((maxY - minY) / 2));

            //assume parabolic dispersion
            var b = maxMin.Max;
            var a = (maxMin.Min - b) / Math.Pow(RadiusAroundShape, 2);


            HashSet<Tuple<int, int>> hashset = new HashSet<Tuple<int, int>>();
            for (int ii = 0; ii < projectedPoints.Count; ii++)
            {
                //prevent adding same points
                if (hashset.Contains(new Tuple<int, int>((int)projectedPoints[ii].X, (int)projectedPoints[ii].Y)))
                {
                    continue;
                }
                else
                {
                    area++;
                    hashset.Add(new Tuple<int, int>((int)projectedPoints[ii].X, (int)projectedPoints[ii].Y));
                    double dist = Math.Sqrt((projectedPoints[ii].X - xcg) * (projectedPoints[ii].X - xcg) +
                                         (projectedPoints[ii].Y - ycg) * (projectedPoints[ii].Y - ycg));

                    double intensity = a * Math.Pow(dist, 2) + b;

                    if ((double)canvas[(int)projectedPoints[ii].X, (int)projectedPoints[ii].Y] + intensity > ushort.MaxValue)
                    {
                        canvas[(int)projectedPoints[ii].X, (int)projectedPoints[ii].Y] = ushort.MaxValue - 1;
                    }
                    else
                    {
                        //AddPSF(canvas, projectedPoints[ii], intensity);
                        canvas[(int)projectedPoints[ii].X, (int)projectedPoints[ii].Y] += (ushort)intensity;
                    }

                }
            }
        }
        private List<PointOnGrid> GetPointsByIndex(List<PointOnGrid> projectedPoints, List<int> newCHull)
        {
            List<PointOnGrid> result = new List<PointOnGrid>();
            for (int ii = 0; ii < newCHull.Count; ii++)
            {
                result.Add(projectedPoints[newCHull[ii]]);
            }

            return result;
        }

        private void AddPSF(ushort[,] canvas, PointOnGrid projectedPoint, double intensity)
        {
            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);


            for (int ii = -1; ii <= 1; ii += 2)
            {
                for (int jj = -1; jj <= 1; jj += 2)
                {
                    int newX = (int)projectedPoint.X + jj;
                    int newY = (int)projectedPoint.Y + ii;

                    bool b1 = newX > 0;
                    bool b2 = newX < width;
                    bool b3 = newY > 0;
                    bool b4 = newY < height;

                    if (b1 && b2 && b3 && b4)
                    {
                        canvas[newX, newY] += (ushort)(intensity / 8); //contribution from eight neighbours
                    }



                }
            }

            canvas[(int)projectedPoint.X, (int)projectedPoint.Y] += (ushort)intensity;

        }

        private List<PointOnGrid> MergeConvexHulls(List<List<PointOnGrid>> convexHulls)
        {
            List<PointOnGrid> result = new List<PointOnGrid>();
            List<PointOnGrid> dbg = new List<PointOnGrid>();


            for (int ii = 0; ii < convexHulls.Count; ii++)
            {
                result.AddRange(convexHulls[ii]);
            }
            if (IsWriteDebug)
                File.WriteAllLines(@"C:\Users\88felix\source\repos\MakeImagesForDescrimination\TempDebug\pointsBefore.csv", result.ToCSV());
            for (int ii = result.Count - 1; ii >= 0; ii--)
            {
                for (int jj = 0; jj < convexHulls.Count; jj++)
                {
                    if (IsPointInPolygon(convexHulls[jj], result[ii]))
                    {
                        if (IsWriteDebug)
                            dbg.Add(result[ii]);
                        result.RemoveAt(ii);
                        break;
                    }
                }


            }

            if (IsWriteDebug)
            {
                File.WriteAllLines(@"C:\Users\88felix\source\repos\MakeImagesForDescrimination\TempDebug\pointsRemoved.csv", dbg.ToCSV());
                File.WriteAllLines(@"C:\Users\88felix\source\repos\MakeImagesForDescrimination\TempDebug\pointsAfter.csv", result.ToCSV());
            }


            return result;
        }

        private void PlotPoints(List<PointOnGrid> partialContour, ref ushort[,] tempArray)
        {
            for (int ii = 0; ii < partialContour.Count; ii++)
            {
                tempArray[(int)partialContour[ii].X, (int)partialContour[ii].Y] = 32000;
            }

            const int icColorCount = 2;
            BitmapSource source;
            int arrayWidth = tempArray.GetLength(1), arrayHeight = tempArray.GetLength(0);

            int iPaddedWidth = arrayWidth;
            if (arrayWidth % 4 != 0)
            {
                iPaddedWidth = (arrayWidth / 4 + 1) * 4;
            }

            Bitmap bmp = new Bitmap(arrayWidth, arrayHeight, PixelFormat.Format16bppGrayScale);
            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadWrite, bmp.PixelFormat);
            int bytesPerPixel = Image.GetPixelFormatSize(bitmapData.PixelFormat) / 8;
            int size = bitmapData.Stride * bitmapData.Height;
            byte[] DataForBitmapOutBytes = new byte[size];
            int col, row, kk;
            try
            {
                // stride is padded to multiple of 4. iterate over the source not destination, otherwise out of range exception

                for (int ii = 0; ii < tempArray.Length; ii++)
                {
                    col = ii % arrayWidth;
                    row = ii / arrayWidth;
                    if (row > 0)
                    {

                    }
                    kk = row * iPaddedWidth * icColorCount + col * icColorCount;
                    DataForBitmapOutBytes[kk] = (byte)(tempArray[row, col]);
                    DataForBitmapOutBytes[kk + 1] = (byte)((tempArray[row, col] & 0xff00) >> 8);


                }
                Marshal.Copy(DataForBitmapOutBytes, 0, bitmapData.Scan0, DataForBitmapOutBytes.Length);
            }
            catch (Exception e)
            {

                throw;
            }
            finally
            {

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

            using (FileStream fs = new FileStream("C:\\Users\\88felix\\Temp\\debuImg" + new Random().Next() + ".png", FileMode.Create, FileAccess.Write))
            {
                PngBitmapEncoder encoder = new PngBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create(source));
                encoder.Save(fs);
            }



        }

        public override uint GetShapeID()
        {
            return ShapeID;
        }
        private List<PointOnGrid> ComputeProjection(StateVector stateVec, List<PointOnGrid> origPoints)
        {   /*
            rotation_matrix = yaw_rotation_matrix * pitch_rotation_matrix * roll_rotation_matrix;
            world_matrix = translation_matrix * rotation_matrix;
            */

            List<PointOnGrid> projPoints = new List<PointOnGrid>();


            double[,] T_mat = Transltions.CreateTranslationMatrix(new double[]
                { stateVec.XCG, stateVec.YCG, stateVec.ZCG });
            double[,] R_Gen = Rotations.CreateRotationThetaPsiPhi(Rotations.ToRadians(stateVec.PitchAngle),
                Rotations.ToRadians(stateVec.YawAngle), Rotations.ToRadians(stateVec.RollAngle));


            double[,] World_mat = BLAS.Multiply(T_mat, R_Gen);

            /*
              [u,v]= Cam_Project_mat*world_matrix*[x y z]
            */

            double[,] Cam_Project_mat = Camera.CreateProjectionMatrix();
            double[,] Launch_to_Camera = Camera.CreateAxisTranformMatrix();
            double[,] FinalMat = BLAS.Multiply(Cam_Project_mat, BLAS.Multiply(Launch_to_Camera, World_mat));


            for (int ii = 0; ii < origPoints.Count; ii++)
            {
                // Apply transformation to original position
                double[] vec = new double[] { origPoints[ii].X, origPoints[ii].Y, origPoints[ii].Z, 1 };
                double[] transf_camera = BLAS.Multiply(FinalMat, vec);


                transf_camera[0] /= transf_camera[2];
                transf_camera[1] /= transf_camera[2];
                transf_camera[2] /= transf_camera[2];

                //double[] transfScaled = BLAS.Multiply(Scales.ScaleMat, transf_camera);
                if (transf_camera[0] > 0 && transf_camera[1] > 0 && transf_camera[0] < Camera.CameraSettings.SensorWidth && transf_camera[1] < Camera.CameraSettings.SensorHeight)
                {
                    projPoints.Add(new PointOnGrid((float)transf_camera[0], (float)transf_camera[1], 0, ii, origPoints[ii].Intensity, origPoints[ii].Component));
                }
            }

            return projPoints;

        }

        public override void Shift(float dx, float dy)
        {
            //double dxWorld = ((Logic.Settings.SensorSize.Width / 2.0 - dx) / Logic.Settings.SensorSize.Width) * stateVec.XCG / 2;
            //double dyWorld = ((Logic.Settings.SensorSize.Height / 2.0 - dy) / Logic.Settings.SensorSize.Height) * stateVec.XCG / 2;
            //stateVec.YCG += dxWorld;
            //stateVec.ZCG += dyWorld;

            stateVec.YCG += dx;
            stateVec.ZCG += dy;
        }
        public override void Rotate(double ang1, double ang2)
        {
            stateVec.PitchAngle = ang1;
            stateVec.YawAngle = ang2;
        }

        public override void Scale(double factor)
        {
            //factor is between 1.5 and 0.5
            // use distance to make object bigger/smaller
            stateVec.XCG *= factor;

            //PointOnGrid cm = new PointOnGrid(0, 0, 0);
            //cm.X = Vertices.Average(pgrid => pgrid.X);
            //cm.Y = Vertices.Average(pgrid => pgrid.Y);
            //cm.Z = Vertices.Average(pgrid => pgrid.Z);
            //for (int ii = 0; ii < Vertices.Count; ii++)
            //{
            //    Vertices[ii] = (Vertices[ii]-cm) * (1 + factor) + cm;
            //}
            //PointOnGrid cm2 = new PointOnGrid(0, 0, 0);
            //cm2.X = Vertices.Average(pgrid => pgrid.X);
            //cm2.Y = Vertices.Average(pgrid => pgrid.Y);
            //cm2.Z = Vertices.Average(pgrid => pgrid.Z);
        }
    

    }
    abstract class baseShape
    {
        public float area { get; set; }

        public List<PointOnGrid> Vertices;
        public abstract uint GetShapeID();
        public double DistanceToClosestContourPoint(Point pnt)
        {
            double minDistance = 1e6;
            for (int ii = 0; ii < Vertices.Count; ii++)
            {
                if (minDistance > Math.Sqrt(Math.Pow(((double)(Vertices[ii].X - pnt.X)), 2) +
                                            Math.Pow(((double)(Vertices[ii].Y - pnt.Y)), 2)))
                {
                    minDistance = Math.Sqrt(Math.Pow(((double)(Vertices[ii].X - pnt.X)), 2) + Math.Pow(((double)(Vertices[ii].Y - pnt.Y)), 2));
                }
            }

            double maxlength = 0;
            for (int ii = 0; ii < Vertices.Count; ii++)
            {
                for (int jj = 0; jj < Vertices.Count; jj++)
                {
                    double distance = Math.Sqrt(Math.Pow(Vertices[ii].X - Vertices[jj].X, 2) +
                                                Math.Pow(Vertices[ii].Y - Vertices[jj].Y, 2));
                    if (distance > maxlength)
                    {
                        maxlength = distance;
                    }
                }
            }
            return maxlength / (minDistance + 0.001);
        }
        public virtual void Shift(float dx, float dy)
        {
            for (int ii = 0; ii < Vertices.Count; ii++)
            {
                Vertices[ii] = new PointOnGrid(Vertices[ii].X + dx, Vertices[ii].Y + dy);
            }
        }

        public virtual void Rotate(double ang1, double ang2)
        {
            Rotate(ang1);
        }
        /// <param name="ang">Degrees of rotation [°]</param>
        public virtual void Rotate(double ang)
        {
            ang = ang * Math.PI / 180;
            double[] RotMat = new double[]
            {
                Math.Cos(ang), -Math.Sin(ang),
                Math.Sin(ang), Math.Cos(ang)
            };

            for (int ii = 0; ii < Vertices.Count; ii++)
            {
                float newX = (float)(Vertices[ii].X * RotMat[0] + Vertices[ii].Y * RotMat[1]);
                float newY = (float)(Vertices[ii].X * RotMat[2] + Vertices[ii].Y * RotMat[3]);
                Vertices[ii] = new PointOnGrid(newX, newY);

            }

            maxXShape = (maxXShape * RotMat[0] + maxYShape * RotMat[1]);
            minXShape = (minXShape * RotMat[0] + minYShape * RotMat[1]);
            maxYShape = (maxXShape * RotMat[2] + maxYShape * RotMat[3]);
            minYShape = (minXShape * RotMat[2] + minYShape * RotMat[3]);
        }
        public virtual void Scale(double factor)
        {
            double avrgX = 0, avrgY = 0;
            int NumOfPoints = Vertices.Count;
            for (int ii = 0; ii < NumOfPoints; ii++)
            {
                avrgX += Vertices[ii].X / NumOfPoints;
                avrgY += Vertices[ii].Y / NumOfPoints;
            }

            for (int ii = 0; ii < NumOfPoints; ii++)
            {
                Vertices[ii] = new PointOnGrid((float)(Vertices[ii].X - avrgX),
                    (float)(Vertices[ii].Y - avrgY));
            }

            for (int ii = 0; ii < NumOfPoints; ii++)
            {
                Vertices[ii] = new PointOnGrid((float)(Vertices[ii].X * factor),
                    (float)(Vertices[ii].Y * factor));
            }

            for (int ii = 0; ii < NumOfPoints; ii++)
            {
                Vertices[ii] = new PointOnGrid((float)(Vertices[ii].X + avrgX),
                    (float)(Vertices[ii].Y + avrgY));
            }

            maxXShape = (maxXShape - avrgX) * factor + avrgX;
            minXShape = (minXShape - avrgX) * factor + avrgX;
            maxYShape = (maxYShape - avrgX) * factor + avrgY;
            minYShape = (minYShape - avrgX) * factor + avrgY;

        }

        internal class MaxMin
        {
            internal ushort Max;
            internal ushort Min;

            public MaxMin(ushort Max, ushort Min)
            {
                this.Max = Max;
                this.Min = Min;
            }
        }

        internal MaxMin oMaxMin;
        internal double maxXShape, maxYShape, minXShape, minYShape;
        public virtual void Render(ref ushort[,] canvas, ref RectangleF rect)
        {
            area = 0;


            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);

            ShiftFromEdge(Vertices, canvas);

            maxXShape = ValidateSize(Vertices.Max(p => p.X), width);
            maxYShape = ValidateSize(Vertices.Max(p => p.Y), height);

            minXShape = ValidateSize(Vertices.Min(p => p.X), width);
            minYShape = ValidateSize(Vertices.Min(p => p.Y), height);

            rect.X = (float)minXShape;
            rect.Y = (float)minYShape;
            rect.Width = (float)(maxXShape - minXShape);
            rect.Height = (float)(maxYShape - minYShape);

            List<double[]> pointsoftarget = new List<double[]>();

            ScaleTargetAndAddToCanvas(canvas, Vertices, oMaxMin, ref rect);

        }

        public void ShiftFromEdge(List<PointOnGrid> vertices, ushort[,] canvas)
        {
            double dx, dy;

            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);

            if (minXShape < 0)
            {
                dx = -minXShape;
            }
            else if (minXShape > width || maxXShape > width)
            {
                dx = width - maxXShape;
            }
            else
            {
                dx = 0;
            }


            if (minYShape < 0)
            {
                dy = -minYShape;
            }
            else if (minYShape > height || maxYShape > height)
            {
                dy = height - maxYShape;
            }
            else
            {
                dy = 0;
            }

            for (int ii = 0; ii < Vertices.Count; ii++)
            {
                Vertices[ii].X += (float)dx;
                Vertices[ii].Y += (float)dy;
            }

            //update maximal positions
            maxXShape += dx;
            minXShape += dx;
            maxYShape += dy;
            minYShape += dy;

        }

        internal double ValidateSize(float val, int max)
        {
            if (val > max)
            {
                val = max;
            }
            else if (val < 0)
            {
                val = 0;
            }
            else
            {
                val = val;
            }

            return val;
        }

        internal virtual void ScaleTargetAndAddToCanvas(ushort[,] canvas, List<PointOnGrid> contour, MaxMin maxMin, ref RectangleF rect)
        {
            int width = canvas.GetLength(1);
            int height = canvas.GetLength(0);



            double xcg = minXShape + (maxXShape - minXShape) / 2;
            double ycg = minYShape + (maxYShape - minYShape) / 2;
            double RadiusAroundShape = Math.Sqrt(((maxXShape - minXShape) / 2) * ((maxXShape - minXShape) / 2) +
                                              ((maxYShape - minYShape) / 2) * ((maxYShape - minYShape) / 2)) ;

            //assume parabolic dispersion
            double b = maxMin.Max;
            double a = (maxMin.Min - b) / Math.Pow(RadiusAroundShape, 2);


            rect.X = (float)minXShape;
            rect.Y = (float)minYShape;
            rect.Width = (float)(maxXShape - minXShape);
            rect.Height = (float)(maxYShape - minYShape);
            
            double dist = 0, intensity = 0;
            int maxComponenets = contour.Max(p => p.Component);
            for (int ii = (int)minYShape; ii < (int)maxYShape; ii++) //Rows
            {
                for (int jj = (int)minXShape; jj < maxXShape; jj++) //Columns
                {
                    PointOnGrid candidate = new PointOnGrid(jj, ii);
                    if (IsPointInPolygon(contour, candidate))
                    {
                        area++;

                        if (true)
                        {
                            dist = DistanceToContour(contour, candidate);
                            double intensity_temp = intensity;
                            intensity = maxMin.Max- dist * (maxMin.Max - maxMin.Min) / (RadiusAroundShape/2) ;
                            var diff = Math.Abs(intensity_temp - intensity);
                        }
                        else if(true)
                        {
                            for (int kk = 1; kk <= maxComponenets; kk++)
                            {
                                if (IsPointInPolygon(contour, candidate))
                                {
                                    //contour.ToCSV().WriteCSV(@"c:\temp\contour.csv");
                                    if (IsPointInPolygon(contour.FindAll(p => p.Component == kk).ToList(), candidate))
                                    {
                                        intensity = maxMin.Min + (maxMin.Max - maxMin.Min) * (1D / maxComponenets) * kk;
                                    }
                                }
                                
                            }

        
                        }
                        else
                        {
                            dist = Math.Sqrt((candidate.X - xcg) * (candidate.X - xcg) +
                                             (candidate.Y - ycg) * (candidate.Y - ycg));
                            intensity = (maxMin.Max - maxMin.Min) - a * Math.Pow(dist, 2);
                        }

                       
                        if ((double)(canvas[(int)candidate.X, (int)candidate.Y] + intensity) > (double)ushort.MaxValue)
                        {
                            canvas[(int)candidate.X, (int)candidate.Y] = ushort.MaxValue - 1;
                        }
                        else
                        {
                            //AddPSF(canvas, projectedPoints[ii], intensity);
                            checked
                            {
                                canvas[(int)candidate.X, (int)candidate.Y] += (ushort)intensity;
                            }

                        }

                    }
                }
            }

        }

        internal bool IsPointInPolygon(List<PointF> polygons, PointF testPoint)
        {
            bool result = false;
            int j = polygons.Count() - 1;
            for (int i = 0; i < polygons.Count(); i++)
            {
                if (polygons[i].Y < testPoint.Y && polygons[j].Y >= testPoint.Y || polygons[j].Y < testPoint.Y && polygons[i].Y >= testPoint.Y)
                {
                    if (polygons[i].X + (testPoint.Y - polygons[i].Y) / (polygons[j].Y - polygons[i].Y) * (polygons[j].X - polygons[i].X) < testPoint.X)
                    {
                        result = !result;
                    }
                }
                j = i;
            }
            return result;
        }
        internal bool IsPointInPolygon(List<PointOnGrid> polygon, PointOnGrid testPoint)
        {
            bool result = false;
            int j = polygon.Count() - 1;
            for (int i = 0; i < polygon.Count(); i++)
            {
                if (polygon[i].Y < testPoint.Y && polygon[j].Y >= testPoint.Y || polygon[j].Y < testPoint.Y && polygon[i].Y >= testPoint.Y)
                {
                    if (polygon[i].X + (testPoint.Y - polygon[i].Y) / (polygon[j].Y - polygon[i].Y) * (polygon[j].X - polygon[i].X) < testPoint.X)
                    {
                        result = !result;
                    }
                }
                j = i;
            }
            return result;
        }

        internal double DistanceToContour(List<PointOnGrid> contourPoints, PointOnGrid point)
        {
            int iDummy = 0;
            return DistanceToContour(contourPoints, point, ref iDummy);
        }
        

        internal double DistanceToContour(List<PointOnGrid> contourPoints, PointOnGrid point,ref int iClosest)
        {
            double minDistance = 1e6;
            
            for (int ii = 0; ii < contourPoints.Count; ii++)
            {
                if (ii == contourPoints.Count - 1)
                {
                    if (contourPoints[ii].DistanceToOtherPoint(contourPoints[0]) > 1e-4)
                    {
                        //todo add intensity effect
                        double dist = GetDistanceToLine(contourPoints[ii], contourPoints[0], point);
                        if (dist < minDistance)
                        {
                            iClosest = ii;
                            minDistance = dist;
                        }
                    }
                    else
                    {
                    }
                }
                else if (ii == 0)
                {
                    if (contourPoints[ii].DistanceToOtherPoint(contourPoints[0]) > 1e-4)
                    {
                        //todo add intensity effect
                        double dist = GetDistanceToLine(contourPoints[ii], contourPoints[0], point);
                        if (dist < minDistance)
                        {
                            iClosest = ii;
                            minDistance = dist;
                        }
                    }
                    else
                    {
                    }
                }
                else
                {
                    if (contourPoints[ii].DistanceToOtherPoint(contourPoints[ii + 1]) > 1e-4)
                    {
                        double dist = GetDistanceToLine(contourPoints[ii], contourPoints[ii + 1], point);

                        if (dist < minDistance)
                        {
                            iClosest = ii;
                            minDistance = dist;
                        }
                    }
                    else
                    {
                    }
                }
            }
            return minDistance;
        }

        private double GetDistanceToLine(PointOnGrid P1, PointOnGrid P2, PointOnGrid point)
        {
            double result = Math.Abs((P2.Y - P1.Y) * point.X - (P2.X - P1.X) * point.Y + P2.X * P1.Y - P2.Y * P1.X) /
                            P1.DistanceToOtherPoint(P2);
            return result;
        }
    }
    [StructLayout(LayoutKind.Explicit)]
    public struct strPoint
    {
        [FieldOffset(0)] public double x;
        [FieldOffset(8)] public double y;
        public strPoint(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
    }
    public class PointOnGrid
    {
        public float Intensity { get; set; }
        public int Index { get; set; }
        public int Component { get; set; }
        public float X { get; set; }
        public float Y { get; set; }
        public float Z { get; set; }
        /// <summary>Initializes a new instance of the class with the specified coordinates. and intensity</summary>
        /// <param name="x">The horizontal position of the point. </param>
        /// <param name="y">The vertical position of the point. </param>
        public PointOnGrid(float x, float y, float z, int index, float intensity, int component)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.Intensity = intensity;
            this.Index = index;
            this.Component = component;
        }
        public PointOnGrid(float x, float y, float z, int index, float intensity)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.Intensity = intensity;
            this.Index = index;
        }
        public PointOnGrid(float x, float y, float z, int index)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.Index = index;
        }
        public PointOnGrid(float x, float y, float z)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
        }
        public PointOnGrid(float x, float y)
        {
            this.X = x;
            this.Y = y;
        }
        public double DistanceToOtherPoint(PointOnGrid otherPoint)
        {
            double result = Math.Sqrt(Math.Pow(Y - otherPoint.Y, 2) - Math.Pow(otherPoint.X - otherPoint.X, 2));
            return result;
        }

        public static PointOnGrid operator /(PointOnGrid p, double scale)
        {
            PointOnGrid res = new PointOnGrid(0, 0, 0, p.Index, p.Intensity);
            res.X = (float)(p.X / scale);
            res.Y = (float)(p.Y / scale);
            res.X = (float)(p.Z / scale);
            return res;
        }
        public static PointOnGrid operator *(PointOnGrid p, double scale)
        {
            PointOnGrid res = new PointOnGrid(0, 0, 0, p.Index, p.Intensity);
            res.X = (float)(p.X * scale);
            res.Y = (float)(p.Y * scale);
            res.Z = (float)(p.Z * scale);
            return res;
        }
        public static PointOnGrid operator +(PointOnGrid p1, PointOnGrid p2)
        {
            PointOnGrid res = new PointOnGrid(0, 0, 0);
            res.X = (float)(p1.X + p2.X);
            res.Y = (float)(p1.Y + p2.Y);
            res.X = (float)(p1.Z + p2.Z);
            return res;
        }
        public static PointOnGrid operator -(PointOnGrid p1, PointOnGrid p2)
        {
            PointOnGrid res = new PointOnGrid(0, 0, 0);
            res.X = (float)(p1.X - p2.X);
            res.Y = (float)(p1.Y - p2.Y);
            res.X = (float)(p1.Z - p2.Z);
            return res;
        }
    }
}