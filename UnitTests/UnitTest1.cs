using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MakeImagesForDescrimination;
using MakeImagesForDiscrimination;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace UnitTests
{
    [TestClass]
    public class UnitTestUtils
    {
        [TestMethod]
        public void TestMethod()
        {
            List<PointOnGrid> lstPoints = new List<PointOnGrid>();
            #region Data
            lstPoints.Add(new PointOnGrid(0.0572F, 1.5740F, 0, 1));
            lstPoints.Add(new PointOnGrid(0.7491F, 1.4529F, 0, 2));
            lstPoints.Add(new PointOnGrid(1.2334F, 1.2972F, 0, 3));
            lstPoints.Add(new PointOnGrid(1.9771F, 1.1761F, 0, 4));
            lstPoints.Add(new PointOnGrid(2.4614F, 1.0724F, 0, 5));
            lstPoints.Add(new PointOnGrid(3.2225F, 0.8994F, 0, 6));
            lstPoints.Add(new PointOnGrid(3.8797F, 0.8129F, 0, 7));
            lstPoints.Add(new PointOnGrid(4.4505F, 1.2626F, 0, 8));
            lstPoints.Add(new PointOnGrid(4.8483F, 1.5567F, 0, 9));
            lstPoints.Add(new PointOnGrid(5.4883F, 2.0928F, 0, 10));
            lstPoints.Add(new PointOnGrid(5.6266F, 0.7783F, 0, 11));
            lstPoints.Add(new PointOnGrid(5.9796F, 0.7846F, 0, 12));
            lstPoints.Add(new PointOnGrid(6.3531F, 0.7783F, 0, 13));
            lstPoints.Add(new PointOnGrid(6.3440F, 1.1836F, 0, 14));
            lstPoints.Add(new PointOnGrid(6.3358F, 1.5394F, 0, 15));
            lstPoints.Add(new PointOnGrid(6.6298F, 1.3145F, 0, 16));
            lstPoints.Add(new PointOnGrid(7.0103F, 1.0897F, 0, 17));
            lstPoints.Add(new PointOnGrid(7.4600F, 0.8821F, 0, 18));
            lstPoints.Add(new PointOnGrid(7.4773F, 0.4324F, 0, 19));
            lstPoints.Add(new PointOnGrid(7.4600F, 0.0692F, 0, 20));
            lstPoints.Add(new PointOnGrid(7.4773F, -0.3632F, 0, 21));
            lstPoints.Add(new PointOnGrid(7.5119F, -0.8302F, 0, 22));
            lstPoints.Add(new PointOnGrid(6.9066F, -1.1416F, 0, 23));
            lstPoints.Add(new PointOnGrid(6.6817F, -0.7610F, 0, 24));
            lstPoints.Add(new PointOnGrid(6.3531F, -0.8129F, 0, 25));
            lstPoints.Add(new PointOnGrid(5.9449F, -0.7597F, 0, 26));
            lstPoints.Add(new PointOnGrid(5.5402F, -0.7783F, 0, 27));
            lstPoints.Add(new PointOnGrid(5.4883F, -1.7296F, 0, 28));
            lstPoints.Add(new PointOnGrid(5.1942F, -1.6431F, 0, 29));
            lstPoints.Add(new PointOnGrid(4.8345F, -1.4191F, 0, 30));
            lstPoints.Add(new PointOnGrid(4.3467F, -1.1243F, 0, 31));
            lstPoints.Add(new PointOnGrid(3.8451F, -0.8475F, 0, 32));
            lstPoints.Add(new PointOnGrid(3.2728F, -0.8812F, 0, 33));
            lstPoints.Add(new PointOnGrid(2.7209F, -1.0378F, 0, 34));
            lstPoints.Add(new PointOnGrid(2.3012F, -1.1067F, 0, 35));
            lstPoints.Add(new PointOnGrid(1.7523F, -1.2453F, 0, 36));
            lstPoints.Add(new PointOnGrid(1.2948F, -1.3497F, 0, 37));
            lstPoints.Add(new PointOnGrid(0.6799F, -1.4356F, 0, 38));
            lstPoints.Add(new PointOnGrid(0.0119F, -1.6258F, 0, 39));
            lstPoints.Add(new PointOnGrid(0.1000F, -1.4000F, 0, 40));
            lstPoints.Add(new PointOnGrid(0.0281F, -1.0200F, 0, 41));
            lstPoints.Add(new PointOnGrid(0.0628F, -0.0136F, 0, 42));
            lstPoints.Add(new PointOnGrid(0.0628F, 0.8366F, 0, 43));


            #endregion Data

            var expectedResult = new[] {39, 1 ,2 ,8, 10 ,18, 19, 22, 23, 28 ,36, 37, 39}; //from matlab

            List<int> indexesOfHull= ConvexHull.GetLocalConvex(lstPoints, 4);
            Debug.Assert(indexesOfHull.Count == expectedResult.Length);

            for (int ii = 0; ii < indexesOfHull.Count; ii++)
            {
                Debug.Assert(indexesOfHull[ii] == expectedResult[ii]-1);// -1 because of Matlab indexing
            }


        }
    }
}

