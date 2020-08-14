using System;
using System.Drawing;
using System.Threading.Tasks;
using System.Reflection;

namespace MakeImagesForDiscrimination
{
    interface IDrawable
    {
        void Render(ref ushort[,] canvas, ref RectangleF rect);
    }

    class DEP : IDrawable
    {
        private void StdCalc(ushort[,] InputArray, ref double Average, ref double StdDev, ref ushort Max, ref ushort Min)
        {
            Average = 0;
            StdDev = 0;
            Max = InputArray[0, 0];
            Min = InputArray[0, 0];
            int ItemsCount = InputArray.Length;
            int ArrayWidth = InputArray.GetLength(1);

            double dSum = 0;
            #region Average/MAX/MIN

            for (int ii = 0; ii < ItemsCount; ii++)
            {

                int row = ii / ArrayWidth;
                int col = ii % ArrayWidth;
                Average += (double)InputArray[row, col] / ItemsCount;
                if (InputArray[row, col] > Max)
                {
                    Max = InputArray[row, col];
                }
                if (InputArray[row, col] < Min)
                {
                    Min = InputArray[row, col];
                }
            }

            #endregion

            #region STDDEV


            for (int ii = 0; ii < ItemsCount; ii++)
            {
                int row = ii / ArrayWidth;
                int col = ii % ArrayWidth;
                double delta = InputArray[row, col] - Average;
                dSum += delta * delta;

            }
            if (ItemsCount > 1)
            {
                StdDev = Math.Sqrt(dSum / (ItemsCount - 1));
            }

            #endregion
        }

        public void Render(ref ushort[,] canvas, ref RectangleF rect)
        {
            Random rand = new Random();
            int x = rand.Next(0, canvas.GetLength(1));
            int y = rand.Next(0, canvas.GetLength(0));
            double avg = 0, std = 0;
            ushort Max = 0, Min = 0;
            StdCalc(canvas, ref avg, ref std, ref Max, ref Min);
            
            //create pixel at random coordinates 4 sigmas below average
            canvas[y, x] = (ushort)((avg - 4 * std) > 0 ? avg - 4 * std : 0);

            x = rand.Next(0, canvas.GetLength(1));
            y = rand.Next(0, canvas.GetLength(0));

            //create pixel at random coordinates 4 sigmas above average
            canvas[y, x] = (ushort)((avg + 4 * std) < UInt16.MaxValue ? avg + 4 * std : UInt16.MaxValue);
        }
    }

    enum eShape
    {
        RV=1,
        Engine=2
    }
    enum eTypeOFImage
    {
        TrainingRV,
        ValidationRV,
        TestRV,
        TrainingEng,
        ValidationEng,
        TestEng
    }

    internal class VersionNumber
    {
        public const string versionNumber = "1.0.2.0";
        public const string versionName = "Make Images for Discrimination: " + versionNumber;
    }

    class Program
    {

        static void Main(string[] args)
        {
            Logic oLogic = new Logic(227, 227);
            oLogic.GenerateImages(200);

        }


    }

    struct LabelData
    {
        public uint image_id;
        public uint category_id;
        public float[] bbox;
        public float area;
        public string originalFile;

        //public float X;
        //public float Y;
        //public float Width;
        //public float Height;
    }
}
