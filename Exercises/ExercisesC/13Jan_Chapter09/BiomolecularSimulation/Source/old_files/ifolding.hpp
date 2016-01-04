/*!*****************************************************************************
 * \file   ifolding.hpp
 * \date   2008-09
 * \author SA and KvM
 * \brief  Defines the interface between the folding simulation and the demo
 *         graphical user interface
 * 
 * Sanne, just derive your simulation class from this interface and implement
 * each of these methods.
 * 
 * class MySim : public ISimInterface
 * {
 * public: 
 *   ... ISimInterface methods ...
 * 
 * protected:
 *   ... your stuff ...
 * };
 ******************************************************************************/

#ifndef __IFOLDING_HPP_
#define __IFOLDING_HPP_

// include our auto-array template class
#include "kvm/autoarray.hpp"

/*!*****************************************************************************
 * \struct TBase
 * \brief  Structure containing position and type of a single protein base
 ******************************************************************************/
struct TBase
{
	int		x;			// x-position on grid
	int		y;			// y-position on grid
	int		t;			// base type
};

// define auto-array for TBase
typedef CAutoArray<TBase> CAABase;

/*!*****************************************************************************
 * \struct TContact
 * \brief  Structure containing ids of connects bases and its associated energy
 ******************************************************************************/
struct TContact
{
	int		  idA;		// first base id (0..nBases-1)
	int		  idB;		// second base id (0..nBases-1)
	double	energy; // contact energy (+1,0,-1)
};

// define auto-array for TContact
typedef CAutoArray<TContact> CAAContact;

// forward declaration of the class
class ISimInterface;

// function to create an instance of this class - provided by the library
ISimInterface* createLatticeSim();

/*!*****************************************************************************
 * \class ISimInterface
 * \brief Interface to the folding simulation
 ******************************************************************************/
class ISimInterface
{
public:

/*!*****************************************************************************
 * \brief Destructor - must be present
 ******************************************************************************/
	virtual
	~ISimInterface	() {}

/*!*****************************************************************************
 * \brief Initialization routine
 ******************************************************************************/
	virtual
	void					init										( ) = 0;

/*!*****************************************************************************
 * \brief De-initialization routine
 ******************************************************************************/
	virtual
	void					free										( ) = 0;

/*!*****************************************************************************
 * \brief Sets the native structure of our model protein
 * \param nBases Number of bases
 * \param pBases Array of data for each base
 ******************************************************************************/
	virtual
	void					setNativeStructure			( CAABase&			rBases ) = 0;

/*!*****************************************************************************
 * \brief Returns the native structure of our model protein
 * \param nBases Number of bases
 * \param pBases Array of data for each base
 ******************************************************************************/
	virtual
	void					getNativeStructure			( CAABase&			rBases ) = 0;

/*!*****************************************************************************
 * \brief Returns the maximum energy for the native structure
 * \retval Maximum energy of the native state
 ******************************************************************************/
	virtual
	double				getMaxEnergy						( ) = 0;

/*!*****************************************************************************
 * \brief Returns the maximum energy for a single bond/contact from the table
 * \retval Maximum energy for a contact
 ******************************************************************************/
	virtual
	double				getMaxContactEnergy			( ) { return 0.0; };

/*!*****************************************************************************
 * \brief Sets the current structure of our model protein
 * \param nBases Number of bases
 * \param pBases Array of data for each base
 ******************************************************************************/
	virtual
	void					setStructure						( CAABase&			rBases ) = 0;

/*!*****************************************************************************
 * \brief Returns the current structure of our model protein
 * \param nBases Number of bases
 * \param pBases Array of data for each base
 * \note  One particle's position is always fixed as a reference point.
 ******************************************************************************/
	virtual
	void					getStructure						( CAABase&			rBases ) = 0;

/*!*****************************************************************************
 * \brief Returns the current structure's contacts
 * \param nContacts Number of contacts
 * \param pContacts Array fo data for each contact
 ******************************************************************************/
	virtual
	void					getStructureAndContacts	( CAABase&			rBases,
																					CAAContact&		rContacts ) = 0;

/*!*****************************************************************************
 * \brief Returns the current state's energy
 * \retval Current state's energy
 * \note The energy depends on both the shape and the code
 ******************************************************************************/
	virtual
	double				getEnergy								( ) = 0;

/*!*****************************************************************************
 * \brief Returns the current state's variation
 * \retval Current state's variation
 * \note The variation depends on the code only
 ******************************************************************************/
	virtual
	double				getVariation						( ) = 0;

/*!*****************************************************************************
 * \brief Returns the current state's degree of folding
 * \retval Current state's degree of folding (0:unfolded, 1:folded)
 ******************************************************************************/
	virtual
	double				getDegreeOfFolding			( ) = 0;

/*!*****************************************************************************
 * \brief Returns the simulation temperature
 * \retval Simulation temperature
 * \todo Define range
 ******************************************************************************/
	virtual
	double				getTemperature					( ) = 0;

/*!*****************************************************************************
 * \brief Sets the simulation temperature
 * \param dTemperature Simulation temperature
 * \todo Define range
 ******************************************************************************/
	virtual
	void					setTemperature					( double				dTemperature ) = 0;

/*!*****************************************************************************
 * \brief Optimizes the code for a given shape
 * \note This should generate a code which is most likely to make the model
 *       protein fold into its native state.
 ******************************************************************************/
	virtual
	void					optimizeCode						( ) = 0;

/*!*****************************************************************************
 * \brief Simulates thermal fluctuations
 * \param nMcSteps Number of MC steps to simulate.
 * \note This simulates the folding process.
 ******************************************************************************/
	virtual
	void					simulate								( int						nMcSteps ) = 0;

/*!*****************************************************************************
 * \brief Stretches the model protein into a horizontal line
 ******************************************************************************/
	virtual
	void					stretch									( ) = 0;

/*!*****************************************************************************
 * \brief Folds the model protein into its native state
 ******************************************************************************/
	virtual
	void					fold										( ) = 0;
};

#endif // __IFOLDING_HPP_
