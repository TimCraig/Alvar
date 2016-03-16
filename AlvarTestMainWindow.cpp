/*
 * AlvarTestMainWindow.cpp
 *
 *  Created on: May 25, 2014
 *      Author: Tim
 */

/*****************************************************************************
 ******************************  I N C L U D E  ******************************
 ****************************************************************************/

#include "AlvarTestMainWindow.h"
#include "DQCVImageUtils.h"
#include "CameraCalibration.h"

#include <QMessageBox>
#include <QStatusBar>
#include <QPainter>
#include <QPen>
#include <QColor>
#include <QFileDialog>
#include <QTextStream>

//#include <Eigen/LU>

/*****************************************************************************
 ***  class AlvarTestMainWindow
 ****************************************************************************/

/*****************************************************************************
 *
 ***  AlvarTestMainWindow::AlvarTestMainWindow
 *
 ****************************************************************************/

AlvarTestMainWindow::AlvarTestMainWindow(
      QWidget* pParent /* = nullptr */) : Base(pParent),
      m_fTagSize(1.0f),
      m_pDetections(nullptr),
      m_pProcessTime(nullptr),
      m_pActionLoadCameraMatrix(nullptr),
      m_strTagSize("TagSize"),
      m_strCameraMatrix("CameraMatrix")
   {
   m_strAppName = "AlvarTest";

   m_ImageSize = QSize(640, 480);
   setGeometry(200, 200, 1400, 800);

   // Just in case it never gets read in
   m_bCameraMatrixValid = false;
   m_CameraMatrix = cv::Mat::eye(3, 3, CV_64F);

   // Until a better way of inputing this is developed
   m_fTagSize = 0.161;

   return;

   } // end of AlvarTestMainWindow::AlvarTestMainWindow

/*****************************************************************************
 *
 ***  AlvarTestMainWindow::Initialize
 *
 ****************************************************************************/

void AlvarTestMainWindow::Initialize(
      Qt::Orientation eOrientation /* = Qt::Horizontal */)
   {
   Base::Initialize(eOrientation);

   SetupStatusBar();

   // Make the captured image valid until one is actually captured
   m_CapturedImage = DCVImage(m_ImageSize.width(), m_ImageSize.height(), CV_8UC3);
   m_CapturedImage.Clear();

   m_Camera.SetRes(m_ImageSize.width(), m_ImageSize.height());
   m_Detector.SetMarkerSize(10);

   m_pCameraHandler->StopCamera();  // TEST!!!
   m_strInputImage = "C:/Projects/Libraries/Alvar/test/benchmark_image2_medium.jpg";
   m_CapturedImage.ReadImage(m_strInputImage.toStdString(), DCVImage::eReadUnchanged);
   ProcessImage(m_CapturedImage);

   return;

   } // end of method AlvarTestMainWindow::Initialize

/*****************************************************************************
 *
 ***  AlvarTestMainWindow::SetupCentralWidget
 *
 ****************************************************************************/

void AlvarTestMainWindow::SetupCentralWidget()
   {
   Base::SetupCentralWidget();

   AppendFileMenu();

   return;

   } // end of method AlvarTestMainWindow::SetupCentralWidget

/******************************************************************************
*
***  AlvarTestMainWindow::AppendFileMenu
*
******************************************************************************/

void AlvarTestMainWindow::AppendFileMenu()
   {
   m_pFileMenu->insertSeparator(m_pActionExit);

   m_pActionLoadCameraMatrix = new QAction(tr("Load Camera Matrix"), this);
   m_pActionLoadCameraMatrix->setStatusTip(
         tr("Load Camera Matrix from Calibration File"));
   m_pFileMenu->insertAction(m_pActionExit, m_pActionLoadCameraMatrix);

   connect(m_pActionLoadCameraMatrix, SIGNAL(triggered()), this,
        SLOT(LoadCameraMatrix()));

   m_pFileMenu->insertSeparator(m_pActionExit);

   return;

   } // end of method AlvarTestMainWindow::AppendFileMenu

/******************************************************************************
*
***  AlvarTestMainWindow::SetupStatusBar()
*
******************************************************************************/

void AlvarTestMainWindow::SetupStatusBar()
   {
   QStatusBar* pSB = statusBar();

   // Tag search time
   QLabel* pTime = new QLabel(tr("Process Time:"));
   pSB->addPermanentWidget(pTime);

   m_pProcessTime = new QLabel(tr("9999"));
   m_pProcessTime->setMinimumSize(m_pProcessTime->sizeHint());
   m_pProcessTime->setAlignment(Qt::AlignRight);
   pSB->addPermanentWidget(m_pProcessTime);

   QLabel* pMS = new QLabel(tr("ms"));
   pSB->addPermanentWidget(pMS);

   m_pProcessTime->setText("");

   // Number of detected Tags
   QLabel* pDetections = new QLabel(tr("Detections:"));
   pSB->addPermanentWidget(pDetections);

   m_pDetections = new QLabel(tr("999999999999999999"));
   m_pDetections->setMinimumSize(m_pDetections->sizeHint());
   m_pDetections->setAlignment(Qt::AlignLeft);
   pSB->addPermanentWidget(m_pDetections);

   m_pDetections->setText("");

   return;

   } // end of method AlvarTestMainWindow::SetupStatusBar

/******************************************************************************
*
***  AlvarTestMainWindow::LoadCameraMatrix
*
******************************************************************************/

void AlvarTestMainWindow::LoadCameraMatrix()
   {
   QString strFileName = QFileDialog::getOpenFileName(this,
         tr("Load Camera Matrix"), QDir::currentPath(),
         tr("Calibration Settings files (*.xml *.yml)"));

   if (!strFileName.isEmpty())
      {
      cv::Mat CameraMatrix(3, 3, CV_64F);
      if (CameraCalibration::ReadCameraMatrix(strFileName.toStdString(),
            CameraMatrix))
         {
         m_CameraMatrix = CameraMatrix;
         m_bCameraMatrixValid = true;
         } // end if
      else
         {
         QString strMsg(tr("Failed to read Settings "));
         strMsg += strFileName;
         QMessageBox::warning(this, tr("Open Settings"), strMsg);
         } // end else
      } // end if

   return;

   } // end of method AlvarTestMainWindow::LoadCameraMatrix

/******************************************************************************
*
***  AlvarTestMainWindow::CameraStarted
*
******************************************************************************/

void AlvarTestMainWindow::CameraStarted()
   {
   Base::CameraStarted();

   return;

   } // end of method AlvarTestMainWindow::CameraStarted

/******************************************************************************
*
***  AlvarTestMainWindow::CameraStopped
*
******************************************************************************/

void AlvarTestMainWindow::CameraStopped()
   {
   Base::CameraStopped();
\
   return;

   } // end of method AlvarTestMainWindow::CameraStopped

/******************************************************************************
*
***  AlvarTestMainWindow::ReprocessImage
*
******************************************************************************/

void AlvarTestMainWindow::ReprocessImage()
   {
   m_CapturedImage.ReadImage(m_strInputImage.toStdString(), DCVImage::eReadUnchanged); // TEST!!!
   Base::ReprocessImage();

   return;

   } // end of method AlvarTestMainWindow::ReprocessImage

/*****************************************************************************
 *
 ***  AlvarTestMainWindow::ProcessImage
 *
 ****************************************************************************/

void AlvarTestMainWindow::ProcessImage(DCVImage& Image)
   {
   m_CapturedImage = Image;
   m_CapturedImage.CvtColor(m_GrayImage, CV_BGR2GRAY);

   m_Timer.start();

   m_Detector.Detect(m_GrayImage, m_Camera, true, true);

   qint64 nElapsed = m_Timer.elapsed();

   m_pProcessTime->setText(QString().setNum(nElapsed));

   DisplayOutput();

   return;

   } // end of method AlvarTestMainWindow::ProcessImage

/******************************************************************************
*
***  AlvarTestMainWindow::DisplayOutput
*
******************************************************************************/

void AlvarTestMainWindow::DisplayOutput()
   {
   const std::vector<alvar::MarkerData>* pMarkers = m_Detector.getMarkers();

   if (pMarkers->size() > 0)
      {
      for (size_t i = 0 ; i < pMarkers->size() ; i++)
         {
         (*pMarkers)[i].Visualize(m_CapturedImage, m_Camera);
         } // end for
      } // end if

   DQImage OutputImage = cvMatToQImage(m_GrayImage);
   DQImage InputImage = cvMatToQImage(m_CapturedImage);

   m_pInputImageWidget->SetImage(InputImage);
   m_pOutputImageWidget->SetImage(OutputImage);

   // Show the number of detected tags in the status bar
   QString strDetections;
   strDetections.setNum(pMarkers->size());
   m_pDetections->setText(strDetections);

   return;

   } // end of method AlvarTestMainWindow::DisplayOutput

/******************************************************************************
*
***  AlvarTestMainWindow::GetTagPose
*
******************************************************************************/

void AlvarTestMainWindow::GetTagPose() const
   {

   return;

   } // end of method AlvarTestMainWindow::GetTagPose

/******************************************************************************
*
***  AlvarTestMainWindow::LoadParameters
*
* Load the last used tag size and camera matrix from the parameter file on
* startup.
*
******************************************************************************/

bool AlvarTestMainWindow::LoadParameters()
   {
   cv::FileStorage FS(GetParameterPath("yml").toStdString(),
         cv::FileStorage::READ);
   bool bRet = FS.isOpened();
   if (bRet)
      {
      cv::FileNode Node = FS[m_strAppName.toStdString()];
      Node[m_strTagSize] >> m_fTagSize;
      Node[m_strCameraMatrix] >> m_CameraMatrix;
      m_bCameraMatrixValid = true;
      } // end if

   return (bRet);

   } // end of method AlvarTestMainWindow::LoadParameters

/******************************************************************************
*
***  AlvarTestMainWindow::SaveParameters
*
* Save the tag size and current camera matrix to the app parameter file on
* exit.
*
******************************************************************************/

bool AlvarTestMainWindow::SaveParameters()
   {
   cv::FileStorage FS(GetParameterPath("yml").toStdString(),
         cv::FileStorage::WRITE);
   bool bRet = FS.isOpened();
   if (bRet)
      {
      FS << m_strAppName.toStdString() << "{"
         << m_strTagSize << m_fTagSize
         << m_strCameraMatrix << m_CameraMatrix
         << "}";
      } // end if

   return (bRet);

   } // end of method AlvarTestMainWindow::SaveParameters


