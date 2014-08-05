/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
  *
 * Copyright (c) 2014 University of Basel
 * All rights reserved.
 */


#ifndef __itkPenalizingMeanSquaresPointSetToImageMetric_h
#define __itkPenalizingMeanSquaresPointSetToImageMetric_h

#include "itkMeanSquaresPointSetToImageMetric.h"


/**
 * This is a simple extension of the itkMeanSquaresPointSetToImageMetric, which adds a regularization
 * of the type regWeight * ||p||^2, where p denote the parameters of the transformation.
 */
namespace itk
{
template< class TFixedPointSet, class TMovingImage >
class ITK_EXPORT PenalizingMeanSquaresPointSetToImageMetric:
  public MeanSquaresPointSetToImageMetric< TFixedPointSet, TMovingImage >
{
public:

  /** Standard class typedefs. */
  typedef PenalizingMeanSquaresPointSetToImageMetric                      Self;
  typedef MeanSquaresPointSetToImageMetric< TFixedPointSet, TMovingImage > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PenalizingMeanSquaresPointSetToImageMetric, Object);

  /** Types transferred from the base class */
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::MeasureType               MeasureType;
  typedef typename Superclass::DerivativeType            DerivativeType;




  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & derivative) const {
       Superclass::GetDerivative(parameters, derivative);
       for( unsigned int i = 0; i < parameters.size(); i++ )
         {
         derivative[i] += parameters[i] * 2 * m_regularizationParameter;
         }
  }

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const {
       double value = Superclass::GetValue(parameters);
       double regValue = 0;
       for( unsigned int par = 0; par < parameters.Size(); par++ )
         {
           regValue += parameters[par] * parameters[par];
         }
       return value + m_regularizationParameter * regValue;
  }

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType & value, DerivativeType & derivative) const {
      Superclass::GetValueAndDerivative(parameters, value, derivative);
      double regValue = 0;
      for( unsigned int i = 0; i < parameters.size(); i++ )
      {
         regValue += parameters[i] * parameters[i];;
         derivative[i] += parameters[i] * 2 * m_regularizationParameter;
      }
      value += m_regularizationParameter * regValue;

  }

  void SetRegularizationParameter(double regParameter) { m_regularizationParameter = regParameter; }

protected:
  PenalizingMeanSquaresPointSetToImageMetric() : m_regularizationParameter(0) {}
  virtual ~PenalizingMeanSquaresPointSetToImageMetric() {}
private:
  PenalizingMeanSquaresPointSetToImageMetric(const Self &); //purposely not implemented
  void operator=(const Self &);                   //purposely not implemented

  double m_regularizationParameter;

};
} // end namespace itk


#endif
