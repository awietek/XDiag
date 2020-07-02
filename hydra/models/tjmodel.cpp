#include <algorithm>
#include <lila/all.h>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indextable.h>
#include <hydra/indexing/indexhubbard.h>
#include <hydra/utils/range.h>
#include <hydra/utils/bitops.h>

#include "tjmodel.h"
#include "hubbardmodeldetail.h"


namespace hydra { namespace models {
    
    template <class coeff_t>
    TJModel<coeff_t>::TJModel
    (BondList bondlist, Couplings couplings, hilbertspaces::hubbard_qn qn)
      : n_sites_(bondlist.n_sites()),
	qn_(qn)
    {
      using hydra::combinatorics::binomial;

      int nup = qn.n_upspins;
      int ndn = qn.n_downspins;
      assert(nup + ndn <= n_sites_);
      dim_ = binomial(n_sites_, nup) * binomial(n_sites_ - nup, ndn);

      // Currently unused operatuors
      std::vector<std::pair<int, int>> currents_;
      std::vector<coeff_t> current_amplitudes_;
      std::vector<std::pair<int, int>> interactions_;
      std::vector<double> interaction_strengths_;
      double U_;

      hubbardmodeldetail::set_hubbard_terms<coeff_t>
      (bondlist, couplings, hoppings_, hopping_amplitudes_,
       currents_, current_amplitudes_,
       interactions_, interaction_strengths_,
       onsites_, onsite_potentials_, 
       szszs_, szsz_amplitudes_,
       exchanges_, exchange_amplitudes_, U_);
    }

    template <class coeff_t>
    lila::Matrix<coeff_t> TJModel<coeff_t>::matrix(bool ninj_term) const
    {
      using state_t = uint64;
      using hydra::combinatorics::up_hole_to_down;
      using hydra::hilbertspaces::Spinhalf;
      using namespace hydra::utils;


      int nup = qn_.n_upspins;
      int ndn = qn_.n_downspins;
      
      // Try allocating upspin / downspin vectors
      std::vector<state_t> upspins;
      std::vector<state_t> dnspins;
      try
	{
	  upspins.resize(dim_);
	  dnspins.resize(dim_);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate upspin/downspin " 
		    << "vectors for TJModel!" << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // Fill vectors holding upspin/downspin configurations
      int64 idx=0;
      auto hs_upspins = Spinhalf<state_t>(n_sites_, nup);
      auto hs_holes_in_ups = Spinhalf<state_t>(n_sites_ - nup, ndn);
      for (state_t ups : hs_upspins)
	for (state_t holes : hs_holes_in_ups)
	  {
	    state_t dns = up_hole_to_down(ups, holes);
	    upspins[idx] = ups;
	    dnspins[idx] = dns;
	    ++idx;
	  }
      assert(idx == dim_);


      // Define lambda function to get indices
      auto index_of_up_dn = 
	[&upspins, &dnspins](state_t const& ups, state_t const& dns)
	{
	  // Binary search the new indices
	  auto up_bounds = std::equal_range(upspins.begin(), 
	  				    upspins.end(), 
	  				    ups);
	  auto up_begin = std::distance(upspins.begin(),
					up_bounds.first); 
	  auto up_end = std::distance(upspins.begin(),
				      up_bounds.second); 
	  auto it = 
	    std::lower_bound(dnspins.begin() + up_begin, 
			     dnspins.begin() + up_end, dns);
	  return std::distance(dnspins.begin(), it);
	};

      // Try allocating the matrix
      lila::Matrix<coeff_t> H;
      try 
	{
	  H.resize(dim_, dim_);
	  Zeros(H);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate matrix for TJModel!" 
		    << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // SzSz terms
      int szsz_idx = 0;
      for (auto pair : szszs_)
	{
	  int s1 = pair.first;
	  int s2 = pair.second;
	  double jz = szsz_amplitudes_[szsz_idx]*0.25;

	  
	  if (std::abs(jz) > 1e-14)
	    {
	      // printf("sites %d %d %d %f\n", s1, s2, szsz_idx, jz);

	      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = upspins[idx];
		  state_t dns = dnspins[idx];
		  
		  bool up1 = gbit(ups, s1);
		  bool up2 = gbit(ups, s2);
		  bool dn1 = gbit(dns, s1);
		  bool dn2 = gbit(dns, s2);

		  if (ninj_term)
		    {
		      if ((up1 && up2) || (dn1 && dn2))
			H(idx, idx) += 0;
		      else if ((up1 && dn2) || (dn1 && up2))
			H(idx, idx) += -2*jz;
		    }
		  else
		    {
		      if ((up1 && up2) || (dn1 && dn2))
			H(idx, idx) += jz;
		      else if ((up1 && dn2) || (dn1 && up2))
			H(idx, idx) += -jz;
		    }
		} 

	    }  // if (std::abs(jz) > 1e-14)
	  ++szsz_idx;
	}  // for (auto pair : szszs_)
      
      // Exchange terms
      int exchange_idx=0;
      for (auto pair: exchanges_)
	{
	  int s1 = std::min(pair.first, pair.second);
	  int s2 = std::max(pair.first, pair.second);
	  coeff_t jx = exchange_amplitudes_[exchange_idx]*0.5;
	  state_t flipmask = ((state_t)1 << s1) | ((state_t)1 << s2);
	  if (std::abs(jx) > 1e-14)
	    {
	      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = upspins[idx];
		  state_t dns = dnspins[idx];
		  bool up1 = gbit(ups, s1);
		  bool up2 = gbit(ups, s2);
		  bool dn1 = gbit(dns, s1);
		  bool dn2 = gbit(dns, s2);
		
		  if ((up1 && dn2) || (dn1 && up2))
		    {
		      state_t flipped_ups = ups ^ flipmask;
		      state_t flipped_dns = dns ^ flipmask;

		      if (popcnt(flipped_ups) != popcnt(ups))
			{
			  using namespace hilbertspaces;
			  std::cout << PrintSpinhalf(64, ups)
				    << std::endl
				    << popcnt(ups) << std::endl
				    << PrintSpinhalf(64, flipped_ups)
				    << std::endl
				    << popcnt(flipped_ups) << std::endl
				    << std::endl;
			}
		      assert(popcnt(flipped_ups) == popcnt(ups));
		      assert(popcnt(flipped_dns) == popcnt(dns));

		      int64 flipped_idx = index_of_up_dn(flipped_ups,
							 flipped_dns);
		      H(flipped_idx, idx) += jx;
		    }

		}  // loop over spin configurations
	    }  // if (std::abs(jz) > 1e-14)
	  ++exchange_idx;
	}  // for (auto pair: exchanges_)


      // Hoppings
      int hopping_idx=0;
      for (auto pair : hoppings_)
      	{
      	  int s1 = std::min(pair.first, pair.second); 
      	  int s2 = std::max(pair.first, pair.second);
      	  coeff_t t = hopping_amplitudes_[hopping_idx];
      	  state_t flipmask = ((state_t)1 << s1) | ((state_t)1 << s2);
	  state_t firstmask = (state_t)1 << s1;
	  
      	  if (std::abs(t) > 1e-14)
      	    {
      	      for (int64 idx=0; idx<dim_; ++idx)
      		{
      		  state_t ups = upspins[idx];
      		  state_t dns = dnspins[idx];
	      
      		  // upspin hopping
      		  if ((dns & flipmask) == 0)
      		    {
      		      if (((ups & flipmask) != 0) &&
			  ((ups & flipmask) != flipmask))
      			{
      			  state_t flipped_ups = ups ^ flipmask;
			  int64 flipped_idx =
			    index_of_up_dn(flipped_ups, dns);
       			  
			  double fermi_up = 
			    popcnt(gbits(ups ^ dns, s2-s1, s1)) % 2==0
			    ? 1. : -1.;

			  if (ups & firstmask)
			    H(flipped_idx, idx) += t * fermi_up;
			  else
			    H(flipped_idx, idx) -= lila::conj(t) * fermi_up;
      			}
      		    }

      		  // dnspin hopping
      		  if ((ups & flipmask) == 0)
      		    {
      		      if (((dns & flipmask) != 0) &&
			  ((dns & flipmask) != flipmask))
      			{
      			  state_t flipped_dns = dns ^ flipmask;
			  int64 flipped_idx =
			    index_of_up_dn(ups, flipped_dns);

			  double fermi_dn = 
      			    popcnt(gbits(ups ^ dns, s2-s1, s1)) % 2==0
			    ? 1. : -1.;

			  if (dns & firstmask)
			    H(flipped_idx, idx) += t * fermi_dn;
			  else
			    H(flipped_idx, idx) -= lila::conj(t) * fermi_dn;
	

      			}
      		    }
      		}  // loop over spin configurations
      	    }  // if (std::abs(t) > 1e-14)
      	  ++hopping_idx;
      	}  // for (auto pair : hoppings_)
    
      return H;

    }

    template class TJModel<double>;
    template class TJModel<complex>;

    template <class coeff_t>
    lila::Matrix<double> TJModel<coeff_t>::szMatrix(int siteIndex) const
    {
      using state_t = uint32;
      using hydra::combinatorics::up_hole_to_down;
      using hydra::hilbertspaces::Spinhalf;
      using namespace hydra::utils;

      assert(0 <= siteIndex && siteIndex < n_sites_); // Check to make sure site index is in range


      int nup = qn_.n_upspins;
      int ndn = qn_.n_downspins;
      
      // Try allocating upspin / downspin vectors
      std::vector<state_t> upspins;
      std::vector<state_t> dnspins;
      try
	{
	  upspins.resize(dim_);
	  dnspins.resize(dim_);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate upspin/downspin " 
		    << "vectors for TJModel!" << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // Fill vectors holding upspin/downspin configurations
      int64 idx=0;
      auto hs_upspins = Spinhalf<state_t>(n_sites_, nup);
      auto hs_holes_in_ups = Spinhalf<state_t>(n_sites_ - nup, ndn);
      for (state_t ups : hs_upspins)
	for (state_t holes : hs_holes_in_ups)
	  {
	    state_t dns = up_hole_to_down(ups, holes);
	    upspins[idx] = ups;
	    dnspins[idx] = dns;
	    ++idx;
	  }
      assert(idx == dim_);


      // Try allocating the matrix
      lila::Matrix<double> sz;
      try 
	{
	  sz.resize(dim_, dim_);
	  Zeros(sz);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate matrix for TJModel!" 
		    << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}
	
     // Iterate through spin configurations

      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = upspins[idx];
		  state_t dns = dnspins[idx];

		  bool up = gbit(ups, siteIndex);
		  bool dn = gbit(dns, siteIndex);

      if (up)
        sz(idx, idx) = 1;
      else if (dn)
        sz(idx, idx) = -1;
      else
        sz(idx, idx) = 0;

    }

    return sz;

  }

    template <class coeff_t>
    lila::Matrix<double> TJModel<coeff_t>::sPlusMatrix(int siteIndex) const
    {
      using state_t = uint32;
      using hydra::combinatorics::up_hole_to_down;
      using hydra::combinatorics::binomial;
      using hydra::hilbertspaces::Spinhalf;
      using namespace hydra::utils;

      assert(0 <= siteIndex && siteIndex < n_sites_); // Check to make sure site index is in range


      int basenup = qn_.n_upspins;
      int basendn = qn_.n_downspins;
      int targetnup = basenup+1;
      int targetndn = basendn-1;

      // Make sure target Hilbert space is valid (return a zero matrix otherwise)
      
      if (targetnup*(targetnup-n_sites_) <= 0 && targetndn*(targetndn-n_sites_) <= 0 ) { 

        int targetdim =  binomial(n_sites_, targetnup) * binomial(n_sites_ - targetnup, targetndn);
        
        // Try allocating upspin / downspin vectors
        std::vector<state_t> baseupspins;
        std::vector<state_t> basednspins;
        std::vector<state_t> targetupspins;
        std::vector<state_t> targetdnspins;
        try
    {
      baseupspins.resize(dim_);
      basednspins.resize(dim_);
      targetupspins.resize(targetdim);
      targetdnspins.resize(targetdim);
    }
        catch(...)
    {
      std::cerr << "Error: Could not allocate upspin/downspin " 
          << "vectors for TJModel!" << std::endl << std::flush;
      exit(EXIT_FAILURE);
    }

        // Fill vectors holding upspin/downspin configurations
        int64 idx=0;
        auto hs_baseupspins = Spinhalf<state_t>(n_sites_, basenup);
        auto hs_holes_in_baseups = Spinhalf<state_t>(n_sites_ - basenup, basendn);
        for (state_t ups : hs_baseupspins)
    for (state_t holes : hs_holes_in_baseups)
      {
        state_t dns = up_hole_to_down(ups, holes);
        baseupspins[idx] = ups;
        basednspins[idx] = dns;
        ++idx;
      }
        assert(idx == dim_);

        idx=0;
        auto hs_targetupspins = Spinhalf<state_t>(n_sites_, targetnup);
        auto hs_holes_in_targetups = Spinhalf<state_t>(n_sites_ - targetnup, targetndn);
        for (state_t ups : hs_targetupspins)
          for (state_t holes : hs_holes_in_targetups)
          {
            state_t dns = up_hole_to_down(ups, holes);
            targetupspins[idx] = ups;
            targetdnspins[idx] = dns;
            ++idx;
          }
          assert(idx == targetdim);


        // Try allocating the matrix
        lila::Matrix<double> sPlus;
        try 
    {
      sPlus.resize(targetdim, dim_);
      Zeros(sPlus);
    }
        catch(...)
    {
      std::cerr << "Error: Could not allocate matrix for TJModel!" 
          << std::endl << std::flush;
      exit(EXIT_FAILURE);
    }
    
        // Define lambda function to get indices for target space
        auto index_of_up_dn = 
    [&targetupspins, &targetdnspins](state_t const& ups, state_t const& dns)
    {
      // Binary search the new indices
      auto up_bounds = std::equal_range(targetupspins.begin(), 
                  targetupspins.end(), 
                  ups);
      auto up_begin = std::distance(targetupspins.begin(),
            up_bounds.first); 
      auto up_end = std::distance(targetupspins.begin(),
                up_bounds.second); 
      auto it = 
        std::lower_bound(targetdnspins.begin() + up_begin, 
             targetdnspins.begin() + up_end, dns);
      return std::distance(targetdnspins.begin(), it);
    };


       // Iterate through spin configurations
       state_t flipmask = ((state_t)1 << siteIndex);

        for (int64 idx=0; idx<dim_; ++idx)
      {
        state_t ups = baseupspins[idx];
        state_t dns = basednspins[idx];

        bool dn = gbit(dns, siteIndex);

        if (dn) {
          state_t flipped_ups = ups ^ flipmask;
          state_t flipped_dns = dns ^ flipmask;
          int64 target_idx = index_of_up_dn(flipped_ups, flipped_dns);
          sPlus(target_idx, idx) = 1;
        }

      }

      return sPlus;

      } else {
        // Return zero matrix

        lila::Matrix<double> sPlus;
      sPlus.resize(1,1);
      Zeros(sPlus);
        return sPlus;
      }

  }

    template <class coeff_t>
    lila::Matrix<double> TJModel<coeff_t>::sMinusMatrix(int siteIndex) const
    {
      using state_t = uint32;
      using hydra::combinatorics::up_hole_to_down;
      using hydra::combinatorics::binomial;
      using hydra::hilbertspaces::Spinhalf;
      using namespace hydra::utils;

      assert(0 <= siteIndex && siteIndex < n_sites_); // Check to make sure site index is in range


      int basenup = qn_.n_upspins;
      int basendn = qn_.n_downspins;
      int targetnup = basenup-1;
      int targetndn = basendn+1;

      if (targetnup*(targetnup-n_sites_) <= 0 && targetndn*(targetndn-n_sites_) <= 0 ) { 

      int targetdim =  binomial(n_sites_, targetnup) * binomial(n_sites_ - targetnup, targetndn);
      
      // Try allocating upspin / downspin vectors
      std::vector<state_t> baseupspins;
      std::vector<state_t> basednspins;
      std::vector<state_t> targetupspins;
      std::vector<state_t> targetdnspins;
      try
	{
	  baseupspins.resize(dim_);
	  basednspins.resize(dim_);
	  targetupspins.resize(targetdim);
	  targetdnspins.resize(targetdim);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate upspin/downspin " 
		    << "vectors for TJModel!" << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // Fill vectors holding upspin/downspin configurations
      int64 idx=0;
      auto hs_baseupspins = Spinhalf<state_t>(n_sites_, basenup);
      auto hs_holes_in_baseups = Spinhalf<state_t>(n_sites_ - basenup, basendn);
      for (state_t ups : hs_baseupspins)
	for (state_t holes : hs_holes_in_baseups)
	  {
	    state_t dns = up_hole_to_down(ups, holes);
	    baseupspins[idx] = ups;
	    basednspins[idx] = dns;
	    ++idx;
	  }
      assert(idx == dim_);

      idx=0;
      auto hs_targetupspins = Spinhalf<state_t>(n_sites_, targetnup);
      auto hs_holes_in_targetups = Spinhalf<state_t>(n_sites_ - targetnup, targetndn);
      for (state_t ups : hs_targetupspins)
        for (state_t holes : hs_holes_in_targetups)
        {
          state_t dns = up_hole_to_down(ups, holes);
          targetupspins[idx] = ups;
          targetdnspins[idx] = dns;
          ++idx;
        }
        assert(idx == targetdim);


      // Try allocating the matrix
      lila::Matrix<double> sMinus;
      try 
	{
	  sMinus.resize(targetdim, dim_);
	  Zeros(sMinus);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate matrix for TJModel!" 
		    << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}
	
      // Define lambda function to get indices for target space
      auto index_of_up_dn = 
	[&targetupspins, &targetdnspins](state_t const& ups, state_t const& dns)
	{
	  // Binary search the new indices
	  auto up_bounds = std::equal_range(targetupspins.begin(), 
	  				    targetupspins.end(), 
	  				    ups);
	  auto up_begin = std::distance(targetupspins.begin(),
					up_bounds.first); 
	  auto up_end = std::distance(targetupspins.begin(),
				      up_bounds.second); 
	  auto it = 
	    std::lower_bound(targetdnspins.begin() + up_begin, 
			     targetdnspins.begin() + up_end, dns);
	  return std::distance(targetdnspins.begin(), it);
	};


     // Iterate through spin configurations
     state_t flipmask = ((state_t)1 << siteIndex);

      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = baseupspins[idx];
		  state_t dns = basednspins[idx];

		  bool up = gbit(ups, siteIndex);

      if (up) {
        state_t flipped_ups = ups ^ flipmask;
        state_t flipped_dns = dns ^ flipmask;
        int64 target_idx = index_of_up_dn(flipped_ups, flipped_dns);
        sMinus(target_idx, idx) = 1;
      }

    }

    return sMinus;

    } else { 

      lila::Matrix<double> sMinus;
	    sMinus.resize(1,1);
	    Zeros(sMinus);
      return sMinus;
    }

    }

  }
}