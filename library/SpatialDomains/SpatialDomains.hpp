////////////////////////////////////////////////////////////////////////////////
//
//  File:  SpatialDomains.hpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Spatial domains definitions and enumerations.
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_SPATIALDOMAINS_H
#define NEKTAR_SPATIALDOMAINS_SPATIALDOMAINS_H

namespace Nektar
{
namespace SpatialDomains
{
/**
 * @brief Indicates the type of element geometry.
 *
 * This property of the element geometry is used to indicate the
 * necessary storage for the element's geometric information and its
 * corresponding computational complexity. In many cases significant
 * savings in both cases can be made based on this information, in
 * comparison to the most generic case.
 */
enum GeomType
{
    eNoGeomType,    ///< No type defined.
    eRegular,       ///< Geometry is straight-sided with constant
                    ///  geometric factors.
    eDeformed,      ///< Geometry is curved or has non-constant factors.
    eMovingRegular, ///< Currently unused.
};

/**
 * @brief Indicates if the geometric information for an element has
 *        been populated.
 */
enum GeomState
{
    eNotFilled, ///< Geometric information has not been generated.
    ePtsFilled  ///< Geometric information has been generated.
};

/**
 * @brief Principle direction for MMF
 */
enum GeomMMF
{
    eTangentX,        ///< X coordinate direction.
    eTangentY,        ///< Y coordinate direction.
    eTangentXY,       ///< XY direction.
    eTangentZ,        ///< Z coordinate direction.
    eTangentCircular, ///< Circular around the centre of domain.
    ePolar,           ///< Polar Coordinate System.
    eSpherical,       ///< Spherical Coordinate System.
    ePseudospherical,       ///< Pseudospherical Coordinate System.
    eLOCAL,           ///< Align the optimal tangent vector.
    eLOCAL1,           ///< Align the optimal tangent vector.
    eLOCAL2,           ///< Align the optimal tangent vector.
    eLOCAL21,           ///< Align the optimal tangent vector.
    eLOCAL3,           ///< Align the optimal tangent vector.
    eLOCALSphere,     ///< Align the second tangent vector.
    eLOCALEllipsoid,     ///< Align the second tangent vector.
};

/*
 * @brief Session file names associated with tangent principle
 * directions.
 */
const char *const GeomMMFMap[] = {"TangentX",  "TangentY",        "TangentXY",
                                  "TangentZ",  "TangentCircular", "Polar",
                                  "Spherical", "Pseudospherical", "LOCAL", "LOCAL1", 
                                  "LOCAL2", "LOCAL21", "LOCAL3", "LOCALSphere","LOCALEllipsoid"};

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_SPATIALDOMAINS_H
