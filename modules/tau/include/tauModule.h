/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2013 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: rsobie                                                   *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#ifndef TAUMODULE_H
#define TAUMODULE_H

#include <framework/core/Module.h>


namespace Belle2 {
  /**
   * \addtogroup modules
   * @{ tauModule @}
   */

  /**
   * tau
   *
   * tau
   *
   */
  class tauModule : public Module {

  public:

    /**
     * Constructor: Sets the description, the properties and the parameters of the module.
     */
    tauModule();

    /**  */
    virtual ~tauModule();

    /**  */
    virtual void initialize();

    /**  */
    virtual void beginRun();

    /**  */
    virtual void event();

    /**  */
    virtual void endRun();

    /**  */
    virtual void terminate();

   // my modules
   void printTracks();
   void printDecay();
 
  private:

  };
}

#endif /* TAUMODULE_H */
