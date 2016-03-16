/*
 * AlvarTestMainWindow.h
 *
 *  Created on: May 25, 2014
 *      Author: Tim
 */

#ifndef ALVARTESTMAINWINDOW_H_
#define ALVARTESTMAINWINDOW_H_

/*****************************************************************************
******************************  I N C L U D E  ******************************
****************************************************************************/

#include "DCVCameraMainWindow.h"
#include "CVImage.h"
#include "DQImage.h"
#include "DQOpenCV.h"

#include <QAction>
#include <QLabel>
#include <QElapsedTimer>

#include "Camera.h"
#include "MarkerDetector.h"

/*****************************************************************************
*
***  class AlvarTestMainWindow
*
*****************************************************************************/

class AlvarTestMainWindow : public DCVCameraMainWindow
   {
      Q_OBJECT

      using Base = DCVCameraMainWindow;

   public:
       explicit AlvarTestMainWindow(QWidget* pParent  = nullptr);
       AlvarTestMainWindow(const AlvarTestMainWindow& src) = delete;

       ~AlvarTestMainWindow() = default;

       AlvarTestMainWindow& operator=(const AlvarTestMainWindow& rhs);

      virtual void Initialize(Qt::Orientation eOrientation = Qt::Horizontal) override;

   protected:
      alvar::MarkerDetector<alvar::MarkerData> m_Detector;
      alvar::Camera m_Camera;

      DCVImage m_GrayImage;

      // Camera matrix as calculated via OpenCV's calibration routines
      cv::Mat m_CameraMatrix;
      // Size of tag black square in user units (for homography)
      float m_fTagSize;

      QElapsedTimer m_Timer;

      // Status bar indicator values
      QLabel* m_pDetections;
      QLabel* m_pProcessTime;

      QAction* m_pActionLoadCameraMatrix;
      bool m_bCameraMatrixValid;

      std::string m_strTagSize;
      std::string m_strCameraMatrix;

      virtual void SetupCentralWidget() override;
      virtual void AppendFileMenu();
      virtual void SetupStatusBar();
      virtual void ProcessImage(DCVImage& Image) override;
      virtual void DisplayOutput();
      virtual void CameraStarted() override;
      virtual void CameraStopped() override;
      virtual void GetTagPose() const;

   protected slots:
      virtual void LoadCameraMatrix();
      virtual void ReprocessImage() override;

      virtual bool LoadParameters() override;
      virtual bool SaveParameters() override;

   private:

   }; // end of class AlvarTestMainWindow



#endif /* ALVARTESTMAINWINDOW_H_ */
