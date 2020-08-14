using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using MakeImagesForDescrimination;
using Newtonsoft.Json;
using PixelFormat = System.Drawing.Imaging.PixelFormat;

namespace MakeImagesForDiscrimination
{
    public class Logic
    {
        private string jsonFile = String.Empty;
        private int[][] palette = new int[ushort.MaxValue][];
        private ushort[,] EOSArray;
        private int NumberOfTotalImages;
        private List<LabelData> lstLabelsTraining = new List<LabelData>();
        private List<LabelData> lstLabelsTest = new List<LabelData>();
        private List<string> lstBackgroundFiles;

        internal class MeshHolder
        {
            public List<PointOnGrid> meshList;
            public string FileName;
        }
        internal static class Settings
        {
            internal static double TrainValidateRatio = 0.8;
            internal static bool UseConcHull = true;

            internal static class ScaleFactors
            {
                internal static double Scale = 1.2;
                internal static double Offset = 0.8;
                internal static int BackgroundScale = 2; //divide by 2^scale
            }
            internal static string PalettePath = @"..\External\Palette\paletteforGalit.txt";
            internal static Dictionary<eTypeOFImage, string> ImagePaths = new Dictionary<eTypeOFImage, string>()
            {
                { eTypeOFImage.TrainingRV,"\\Training\\RV\\"},
                {eTypeOFImage.ValidationRV,"\\Validation\\RV\\"},
                {eTypeOFImage.TestRV,"\\Test\\RV\\"},
                {eTypeOFImage.TrainingEng,"\\Training\\Eng\\"},
                {eTypeOFImage.ValidationEng,"\\Validation\\Eng\\"},
                {eTypeOFImage.TestEng,"\\Test\\Eng\\"}
            };
            internal static class Geometries
            {
                public static int RVMinH = 20;
                public static int RVMaxH = 40;
                public static int RVMinW = 5;
                public static int RVMaxW = 10;
                public static int EngMinH = 10;
                public static int EngMaxH = 20;
                public static int EngMinW = 45;
                public static int EngMaxW = 50;
                static Random rnd = new Random();
                public static int GetRVRandomH()
                {
                    return rnd.Next(RVMinH, RVMaxH);
                }
                public static int GetRVRandomW()
                {
                    return rnd.Next(RVMinW, RVMaxW);
                }
                public static int GetEngRandomH()
                {
                    return rnd.Next(EngMinH, EngMaxH);
                }
                public static int GetEngRandomW()
                {
                    return rnd.Next(EngMinW, EngMaxW);
                }
            }
            internal static class Temperatures
            {
                internal static baseShape.MaxMin RvTemp = new baseShape.MaxMin(35000, 20000);
                internal static baseShape.MaxMin EngTemp = new baseShape.MaxMin(50000, 20000);
                internal static baseShape.MaxMin Bckgrnd = new baseShape.MaxMin(9100, 8900);
            }

            internal static Size SensorSize;

            static Settings()
            {

            }

            public static void LoadMeshes()
            {
                FillPoints(eTypeOFImage.TestEng, ref TestEngMeshPoints);
                FillPoints(eTypeOFImage.TestRV, ref TestRVMeshPoints);
                FillPoints(eTypeOFImage.ValidationEng, ref ValidateEngMeshPoints);
                FillPoints(eTypeOFImage.ValidationRV, ref ValidateRVMeshPoints);
                FillPoints(eTypeOFImage.TrainingEng, ref TrainEngMeshPoints);
                FillPoints(eTypeOFImage.TrainingRV, ref TrainRVMeshPoints);
            }
            private static void FillPoints(eTypeOFImage typeOfImage, ref MeshHolder[] meshpoints)
            {
                string[] s1MeshFiles;
                s1MeshFiles = Directory.GetFiles(dictPaths[typeOfImage], "*.txt");
                meshpoints = new MeshHolder[s1MeshFiles.Length];
                for (int jj = 0; jj < s1MeshFiles.Length; jj++)
                {
                    meshpoints[jj] = new MeshHolder();
                    meshpoints[jj].meshList = new List<PointOnGrid>();
                    meshpoints[jj].FileName = s1MeshFiles[jj];
                    LoadMeshFile(s1MeshFiles[jj], ref meshpoints[jj].meshList);
                }
            }

            static Dictionary<eTypeOFImage, string> dictPaths = new Dictionary<eTypeOFImage, string>()
            {
                {eTypeOFImage.TestEng,@"..\External\Shapes\Test\Eng"} ,
                {eTypeOFImage.TestRV,@"..\External\Shapes\Test\RV"} ,
                {eTypeOFImage.ValidationEng,@"..\External\Shapes\Validate\Eng"} ,
                {eTypeOFImage.ValidationRV,@"..\External\Shapes\Validate\RV"} ,
                {eTypeOFImage.TrainingEng,@"..\External\Shapes\Train\Eng"} ,
                {eTypeOFImage.TrainingRV,@"..\External\Shapes\Train\RV"} ,

            };

            public static MeshHolder[] TrainRVMeshPoints;
            public static MeshHolder[] TrainEngMeshPoints;
            public static MeshHolder[] TestRVMeshPoints;
            public static MeshHolder[] TestEngMeshPoints;
            public static MeshHolder[] ValidateRVMeshPoints;
            public static MeshHolder[] ValidateEngMeshPoints;

            private static bool LoadMeshFile(string fpath, ref List<PointOnGrid> meshPoints)
            {
                bool result = true;
                string[] s1ShapeVertices = File.ReadAllLines(fpath, Encoding.ASCII);
                meshPoints = new List<PointOnGrid>();
                string[] s1Values;
                float y, z, x;
                int componenet;

                try
                {
                    for (int ii = 0; ii < s1ShapeVertices.Length; ii++)
                    {
                        s1Values = s1ShapeVertices[ii].Split(new char[] { ' ', ',' });
                        y = float.Parse(s1Values[0]);
                        z = float.Parse(s1Values[1]);
                        x = float.Parse(s1Values[2]);
                        componenet = int.Parse(s1Values[3]);
                        meshPoints.Add(new PointOnGrid(x, y, z, ii, 5, componenet));
                    }

                }
                catch (Exception ex)
                {
                    result = false;
                }

                return result;
            }
        }

        private string TimeOfRun;

        #region Init
        public Logic(int Height, int Width)
        {
            Settings.LoadMeshes();
            this.Width = Width;
            this.Height = Height;
            Settings.SensorSize = new Size(Width, Height);
        }
        public Logic()
        {
            Width = 224;
            Height = 224;
            Settings.SensorSize = new Size(Width, Height);
        }
        private int Width;
        private int Height;
        #endregion Init

        public bool GenerateImages()
        {
            return GenerateImages(200);
        }

        public bool GenerateImages(int numberOfImages)
        {
            NumberOfTotalImages = numberOfImages;
            bool bRun = true, bSuccess = false;
            Stopwatch sw = Stopwatch.StartNew();
            //palette = LoadPalette(@"..\External\Palette\hot16Palette.txt");
            palette = Utils.LoadPalette(Settings.PalettePath);

            Random rnd = new Random(DateTime.UtcNow.Millisecond);
            TimeOfRun = DateTime.Now.ToString("HHmmss.fff");
            jsonFile = TimeOfRun + "_" + rnd.Next(0, 666).ToString("000") + ".jsn";



            int GlobalFrmCnt = 0;
            //For NumberOfTotalImages
            lstLabelsTraining.Clear();

            int ii = 0;
            while (ii < Settings.TrainValidateRatio * NumberOfTotalImages)
            {
                //RV creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.TrainingRV);
                // Eng creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.TrainingEng);
                ii++;
            }

            //For Validation
            ii = 0;
            for (ii = 0; ii < (1 - Settings.TrainValidateRatio) * NumberOfTotalImages; ii++)
            {
                //RV creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.ValidationRV);
                // Eng creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.ValidationEng);
                ii++;
            }

            GlobalFrmCnt = 0;
            //For Test
            ii = 0;
            while (ii < (1 - Settings.TrainValidateRatio) * NumberOfTotalImages)
            {
                //RV creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.TestRV);
                // Eng creation
                GenerateShape(ref GlobalFrmCnt, ii, eTypeOFImage.TestEng);
                ii++;
            }

            WriteToJsonFile(lstLabelsTraining, "Train");
            WriteToJsonFile(lstLabelsTest, "Test");
            return bSuccess;
        }

        private void GenerateShape(ref int GlobalFrmCnt, int ii, eTypeOFImage eType)
        {
            string sConsoleMessage = null;
            baseShape shape = null;
            List<LabelData> localLabelsList = null;
            string sOriginaleFileName = "na";
            MeshHolder meshHolder;
            double ratio = 1;
            
            switch (eType)
            {
                case eTypeOFImage.TestEng:
                    sConsoleMessage = "Created Engine for Test image ";
                    //shape = new Cylinder(Settings.Geometries.GetEngRandomH(), Settings.Geometries.GetEngRandomW());
                    meshHolder = Settings.TestEngMeshPoints.RandomIndex();
                    shape = new ComplexShape(meshHolder.meshList, Settings.Temperatures.EngTemp, eShape.Engine);
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTest;
                    ratio = (1 - Settings.TrainValidateRatio);
                    break;
                case eTypeOFImage.TestRV:
                    sConsoleMessage = "Created RV for Test image ";
                    meshHolder = Settings.TestRVMeshPoints.RandomIndex();
                    shape = new ComplexShape(meshHolder.meshList, Settings.Temperatures.RvTemp, eShape.RV);
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTest;
                    ratio = (1 - Settings.TrainValidateRatio);
                    break;
                case eTypeOFImage.TrainingEng:
                    sConsoleMessage = "Created Engine for Training image ";
                    meshHolder = Settings.TrainEngMeshPoints.RandomIndex();
                    shape = new ComplexShape(meshHolder.meshList, Settings.Temperatures.EngTemp, eShape.Engine);
                    //shape = new Cylinder(Settings.Geometries.GetEngRandomH(), Settings.Geometries.GetEngRandomW());
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTraining;
                    ratio = Settings.TrainValidateRatio;
                    break;
                case eTypeOFImage.TrainingRV:
                    sConsoleMessage = "Created RV for Training image ";
                    //shape = new Conus(Settings.Geometries.GetRVRandomH(), Settings.Geometries.GetRVRandomW());
                    meshHolder = Settings.TrainRVMeshPoints.RandomIndex();
                    shape = new ComplexShape(Settings.TrainRVMeshPoints.RandomIndex().meshList, Settings.Temperatures.RvTemp, eShape.RV);
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTraining;
                    ratio = Settings.TrainValidateRatio;
                    break;
                case eTypeOFImage.ValidationEng:
                    sConsoleMessage = "Created Engine for Validation image ";
                    //shape = new Cylinder(Settings.Geometries.GetEngRandomH(), Settings.Geometries.GetEngRandomW());
                    meshHolder = Settings.ValidateEngMeshPoints.RandomIndex();
                    shape = new ComplexShape(meshHolder.meshList, Settings.Temperatures.EngTemp, eShape.Engine);
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTraining;
                    ratio = (1 - Settings.TrainValidateRatio);
                    break;
                case eTypeOFImage.ValidationRV:
                    sConsoleMessage = "Created RV for Validation image ";
                    //shape = new Conus(Settings.Geometries.GetRVRandomH(), Settings.Geometries.GetRVRandomW());
                    meshHolder = Settings.ValidateRVMeshPoints.RandomIndex();
                    shape = new ComplexShape(meshHolder.meshList, Settings.Temperatures.RvTemp, eShape.RV);
                    sOriginaleFileName = meshHolder.FileName;
                    localLabelsList = lstLabelsTraining;
                    ratio = (1 - Settings.TrainValidateRatio);
                    break;

            }

            if (shape != null && sConsoleMessage != null)
            {
                Random rnd = new Random(DateTime.UtcNow.Millisecond);

                shape.Scale(Settings.ScaleFactors.Scale * rnd.NextDouble() + Settings.ScaleFactors.Offset);
                //dirty hack, i rely on the fact that only ComplexShape implements 2 arguments
                shape.Rotate(rnd.Next(1, 359), rnd.Next(1, 359));
                shape.Shift(rnd.Next(-5, 5), rnd.Next(-5, 5));

                RectangleF rect = RectangleF.Empty;
                EOSArray = new ushort[Height, Width];

                //Choose random background (both real and generated)
                if (rnd.Next(0, 3) > 0)
                {
                    AddRealBackGround(ref EOSArray);
                }
                else
                {
                    AddRandomBackGround(ref EOSArray);
                }

                shape.Render(ref EOSArray, ref rect);
                DEP dep = new DEP();
                RectangleF dummyRect = RectangleF.Empty;
                dep.Render(ref EOSArray, ref dummyRect);

                if (rect != RectangleF.Empty)
                {
                    CreateLabel(GlobalFrmCnt, shape.GetShapeID(), rect, shape.area, localLabelsList, sOriginaleFileName);
                    SaveToImage(eType, ref EOSArray, GlobalFrmCnt.ToString("0000"));
                    Console.WriteLine(sConsoleMessage + ii + " out of " + ratio * NumberOfTotalImages);
                    GlobalFrmCnt++;
                }
                else
                {
                    GenerateShape(ref GlobalFrmCnt,ii, eType);
                }
            }
        }

        private void CreateLabel(int GlobalFrmCnt, uint targetID, RectangleF rect, float inputArea, List<LabelData> lstLabels, string sOrigFile)
        {
            LabelData lbl = new LabelData();
            lbl.image_id = (ushort)GlobalFrmCnt;
            lbl.category_id = targetID;
            lbl.bbox = new float[] { rect.X, rect.Y, rect.Width, rect.Height };
            lbl.area = inputArea;
            lbl.originalFile = sOrigFile;
            lstLabels.Add(lbl);
        }


        private void WriteToJsonFile(List<LabelData> lstOfLabels, string suffix)
        {
            string JsonPath = "Results_" + TimeOfRun + "\\";
            if (!Directory.Exists(JsonPath))
            {
                Directory.CreateDirectory(JsonPath);
            }

            string content = JsonConvert.SerializeObject(lstOfLabels, Formatting.Indented);

            FileStream fs;
            StreamWriter streamWriter;
            string jsonFileWithSuffix = jsonFile.Insert(jsonFile.Length - 4, "_" + suffix);
            if (!File.Exists(JsonPath + jsonFileWithSuffix))
            {
                using (fs = new FileStream(JsonPath + jsonFileWithSuffix, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite))
                using (streamWriter = new StreamWriter(fs))
                {
                    streamWriter.WriteLine(content);
                }


            }
            else
            {
                using (fs = new FileStream(JsonPath + jsonFileWithSuffix, FileMode.Append, FileAccess.Write, FileShare.ReadWrite))
                using (streamWriter = new StreamWriter(fs))
                {
                    streamWriter.WriteLine(content);
                }

            }
        }

        private void SaveToImage(eTypeOFImage typeofimage, ref ushort[,] InputEOSArray, string FileNum)
        {
            string localPath = AppDomain.CurrentDomain.BaseDirectory;
            string resultsRawPath = localPath + "Results_" + TimeOfRun + "\\Raw" + Settings.ImagePaths[typeofimage];
            Bitmap bmp;
            BitmapSource bmpSourceGray;
            if (!Directory.Exists(resultsRawPath))
            {
                Directory.CreateDirectory(resultsRawPath);
            }

            string resultsBMPPath = localPath + "\\Results_" + TimeOfRun + "\\Bmp" + Settings.ImagePaths[typeofimage];
            if (!Directory.Exists(resultsBMPPath))
            {
                Directory.CreateDirectory(resultsBMPPath);
            }


            string resultsPNGPath = localPath + "\\Results_" + TimeOfRun + "\\PNG" + Settings.ImagePaths[typeofimage];
            if (!Directory.Exists(resultsPNGPath))
            {
                Directory.CreateDirectory(resultsPNGPath);
            }

            int Width = InputEOSArray.GetLength(1);
            int Height = InputEOSArray.GetLength(0);


            bmpSourceGray = Utils.ArrayToGrayBitmap(EOSArray);
            bmpSourceGray.SaveTo16bPNG(resultsPNGPath + "img_" + FileNum + ".png");

            bmp = Utils.ArrayToBitmap(EOSArray, palette);
            bmp.Save(resultsBMPPath + "img_" + FileNum + ".bmp", ImageFormat.Bmp);

            byte[] by1EOSArray = new byte[EOSArray.Length * sizeof(ushort)];
            Buffer.BlockCopy(EOSArray, 0, by1EOSArray, 0, by1EOSArray.Length);
            File.WriteAllBytes(resultsRawPath + Width + "x" + Height + "_img_" + FileNum + ".raw", by1EOSArray);
        }



        private void AddRandomBackGround(ref ushort[,] InputEosArray)
        {
            int Height = InputEosArray.GetLength(0), Width = InputEosArray.GetLength(1);
            Random rand = new Random();
            ushort[] EOSArrayBackround = new ushort[Height * Width];

            for (int ii = 0; ii < EOSArrayBackround.Length; ii++)
            {
                EOSArrayBackround[ii] =
                    (ushort)rand.Next(Settings.Temperatures.Bckgrnd.Min, Settings.Temperatures.Bckgrnd.Max);
            }
            Buffer.BlockCopy(EOSArrayBackround, 0, EOSArray, 0, EOSArrayBackround.Length * sizeof(ushort));
        }


        private void AddRealBackGround(ref ushort[,] InputEosArray)
        {
            int Height = InputEosArray.GetLength(0), Width = InputEosArray.GetLength(1);
            string localPath = AppDomain.CurrentDomain.BaseDirectory;
            string bckgrndsRawPath = localPath + "..\\External\\RawBackgroundNoise\\";
            if (lstBackgroundFiles == null)
            {
                lstBackgroundFiles = new List<string>(Directory.GetFiles(bckgrndsRawPath, "*.raw"));
            }

            Random rand = new Random();
            byte[] by1FullWindowBckgrnd = File.ReadAllBytes(lstBackgroundFiles[rand.Next(0, lstBackgroundFiles.Count)]);
            ushort[,] us1FullWindowBckgrnd = new ushort[1024, 1024];
            Buffer.BlockCopy(by1FullWindowBckgrnd, 0, us1FullWindowBckgrnd, 0, by1FullWindowBckgrnd.Length);

            int xc = us1FullWindowBckgrnd.GetLength(1) / 2;
            int yc = us1FullWindowBckgrnd.GetLength(0) / 2;
            int xLeft = xc - Width / 2;
            int yTop = yc - Height / 2;
            int offset = ((yTop - 1) * Width + xLeft) * sizeof(ushort); //offset in bytes in original array
            for (int ii = 0; ii < Height; ii++)
            {
                Buffer.BlockCopy(us1FullWindowBckgrnd, offset + (ii * 1024 * sizeof(ushort)), InputEosArray, ii * Width * sizeof(ushort), Width * sizeof(ushort));
            }

            #region ScaleNoise to ~8k adu

            for (int ii = 0; ii < Width * Height; ii++)
            {
                int col = ii % Width;
                int row = ii / Width;

                EOSArray[row, col] = (ushort)(EOSArray[row, col] >> Settings.ScaleFactors.BackgroundScale);
            }

            #endregion ScaleNoise to ~8k adu
        }
    }
}