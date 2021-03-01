// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
using namespace Rcpp;

#include <bitset>
#include <RcppThread.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
int sample_func(const Rcpp::NumericMatrix &r, int ix, int iy, int lookup_value) {
  const int dx = 300;
  
  int count = 0;
  
  
  for (int y=iy-dx; y<iy+dx; ++y)
    for (int x=ix-dx; x<ix+dx; ++x) 
      if ( x>=0 && x<=r.ncol() && y>=0 && y<=r.nrow() && r(y, x) == lookup_value ) 
        ++count;
  return count;
  
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_func2(const Rcpp::NumericMatrix &r, Rcpp::NumericVector xidx, Rcpp::NumericVector yidx, int lookup_value) {
  const int dx = 300;
  
  Rcpp::NumericVector result( xidx.length() );
  
  for (int i = 0; i<xidx.length(); ++i)  {
    int count = 0;
    int ix = xidx(i);
    int iy = yidx(i);

    for (int y=iy-dx; y<iy+dx; ++y)
      for (int x=ix-dx; x<ix+dx; ++x) 
        if ( x>=0 && x<=r.ncol() && y>=0 && y<=r.nrow() && r(y, x) == lookup_value ) 
          ++count;
    
    result(i) = count;
  }
  return result;
      
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_func3(const Rcpp::NumericMatrix &r, 
                                 Rcpp::NumericVector xidx, 
                                 Rcpp::NumericVector yidx, 
                                 Rcpp::NumericVector lookup_value,
                                 Rcpp::NumericVector xdelta,
                                 Rcpp::NumericVector ydelta) {

  Rcpp::NumericVector result( xidx.length() );
  
  for (int i = 0; i<xidx.length(); ++i)  {
    int count = 0;
    int startx = xidx(i);
    int starty = yidx(i);
    int value = lookup_value(i);
    
    for (int j = 0; j< xdelta.length(); ++j) {
      int x = startx + xdelta(j);
      int y = starty + ydelta(j);
      if ( x>=0 && x<=r.ncol() && y>=0 && y<=r.nrow() && r(y, x) == value ) 
        ++count;
    }
    result(i) = count;
  }
    
  return result;
    
  
}



// [[Rcpp::export]]
Rcpp::IntegerVector countPxCentroid(const Rcpp::IntegerMatrix &r, 
                               Rcpp::NumericVector xidx, 
                               Rcpp::NumericVector yidx, 
                               Rcpp::NumericVector lookup_value,
                                 Rcpp::NumericVector xdelta,
                                 Rcpp::NumericVector ydelta) {
  
  Rcpp::IntegerVector result( xidx.length() );
  
  int ncol = r.ncol();
  int nrow = r.nrow();
  int dlength = xdelta.length();
  int n=0;
  
  for (int i = 0; i<xidx.length(); ++i) {
    int value = lookup_value(i);
    int x = xidx(i);
    int y = yidx(i);
    int n=0;
    for (int j = 0; j< dlength; ++j) {
      int dx = x + xdelta(j);
      int dy = y + ydelta(j);
      if ( dx>=0 && dx<ncol && dy>=0 && dy<nrow && r(dy, dx) == value ) 
        n++;
    }
    result(i) = n;
    if (i % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      Rcout << i << "/" << xidx.length() << " patches" << std::endl;
    }
  } 
  Rcout << "Finsihed! " << std::endl;
  
  return result;

}

// [[Rcpp::export]]
Rcpp::IntegerMatrix countPxCentroid2(const Rcpp::IntegerMatrix &r, 
                                     Rcpp::NumericVector xidx, 
                                     Rcpp::NumericVector yidx, 
                                     Rcpp::NumericVector lookup_value,
                                     Rcpp::NumericVector xdelta,
                                     Rcpp::NumericVector ydelta) {
  
  Rcpp::IntegerMatrix result( xidx.length(),4 );
  
  int ncol = r.ncol();
  int nrow = r.nrow();
  int dlength = xdelta.length();
  
  for (int i = 0; i<xidx.length(); ++i) {
    int value = lookup_value(i); // the year we are looking for
    int x = xidx(i); // starting point of patch i
    int y = yidx(i);
    int n=0; // count equal values
    int n_total = 0; // count all disturbance pixel
    int n_prev = 0;
    int n_next = 0;
    for (int j = 0; j< dlength; ++j) {
      int dx = x + xdelta(j);
      int dy = y + ydelta(j);
      if ( dx>=0 && dx<ncol && dy>=0 && dy<nrow ) {
        int pxval = r(dy,dx);
        if (pxval == value)
          ++n;
        if (pxval > 0)
          ++n_total;
        if (pxval == value + 1)
          ++n_prev;
        if (pxval == value - 1)
          ++n_next;
      } 
      
    }
    result(i,0) = n;
    result(i,1) = n_total;
    result(i,2) = n_prev;
    result(i,3) = n_next;
    if (i % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      Rcout << i << "/" << xidx.length() << " patches" << std::endl;
    }
  } 
  Rcout << "Finsihed! " << std::endl;
  
  return result;
  
}


/**
 * Inputs:
 * r: disturbance map (as matrix)
 * xidx/yidx: col/row of the disturbance patch (vector with N patches)
 * lookup_value: value (year) (vector with N patches)
 * xdelta/ydelta: offsets of the kernel
 * 
 * Result: Matrix with 5 columns:
 * (1) count px of the same year
 * (2) total count of disturbed px over the 30years
 * (3) count px of the previous year
 * (4) count px of the next year
 * (5) total number of patches with the same year
 */
// [[Rcpp::export]]
Rcpp::IntegerMatrix countPxCentroid3(const Rcpp::IntegerMatrix &r, 
                                     Rcpp::NumericVector xidx, 
                                     Rcpp::NumericVector yidx, 
                                     Rcpp::NumericVector lookup_value,
                                     Rcpp::NumericVector xdelta,
                                     Rcpp::NumericVector ydelta,
                                     const Rcpp::NumericVector &patchdata) {
  
  Rcpp::IntegerMatrix result( xidx.length(), 5 );
  
  int ncol = r.ncol();
  int nrow = r.nrow();
  int dlength = xdelta.length();
  
  
  std::bitset<100000> patches; // number of bits must be higher than the largest patchId
  
  for (int i = 0; i<xidx.length(); ++i) {
    int value = lookup_value(i); // the year we are looking for
    int x = xidx(i); // starting point of patch i
    int y = yidx(i);
    int n=0; // count equal values
    int n_total = 0; // count all disturbance pixel
    int n_prev = 0;
    int n_next = 0;
    
    patches.reset(); // clear list of patches
    
    for (int j = 0; j< dlength; ++j) {
      int dx = x + xdelta(j);
      int dy = y + ydelta(j);
      if ( dx>=0 && dx<ncol && dy>=0 && dy<nrow ) {
        // pixel values
        int pxval = r(dy,dx);
        if (pxval == value)
          ++n;
        if (pxval > 0)
          ++n_total;
        if (pxval == value - 1)
          ++n_prev;
        if (pxval == value + 1)
          ++n_next;
        // patches: flag each ID with a single bit
        int patch_val = patchdata(dy * ncol + dx);
        if (patch_val > 0 )
          patches.set(patch_val);
        
      } 
      
    }
    result(i,0) = n;
    result(i,1) = n_total;
    result(i,2) = n_prev;
    result(i,3) = n_next;
    result(i,4) = patches.count();
    
    if (i % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      Rcout << i << "/" << xidx.length() << " patches" << std::endl;
    }
  } 
  Rcout << "Finsihed! " << std::endl;
  
  return result;
  
}


/**
 * Inputs:
 * r: disturbance map (as matrix)
 * xidx/yidx: col/row of the disturbance patch (vector with N patches)
 * lookup_value: value (year) (vector with N patches)
 * xdelta/ydelta: offsets of the kernel
 * 
 * Result: Matrix with 5 columns:
 * (1) count px of the same year
 * (2) total count of disturbed px over the 30years
 * (3) count px of the previous year
 * (4) count px of the next year
 * (5) total number of patches with the same year
 */
// [[Rcpp::export]]
Rcpp::IntegerMatrix countPxCentroidMT(const Rcpp::IntegerMatrix &r, 
                                     Rcpp::NumericVector xidx, 
                                     Rcpp::NumericVector yidx, 
                                     Rcpp::NumericVector lookup_value,
                                     Rcpp::NumericVector xdelta,
                                     Rcpp::NumericVector ydelta,
                                     const Rcpp::NumericVector &patchdata) {
  
  Rcpp::IntegerMatrix result( xidx.length(), 5 );
  
  int ncol = r.ncol();
  int nrow = r.nrow();
  int dlength = xdelta.length();
  
  
  std::bitset<1000000> patches; // number of bits must be higher than the largest patchId
  
  //for (int i = 0; i<xidx.length(); ++i) {
  RcppThread::parallelFor(0, xidx.length(), [&](int i) {  
    int value = lookup_value(i); // the year we are looking for
    int x = xidx(i); // starting point of patch i
    int y = yidx(i);
    int n=0; // count equal values
    int n_total = 0; // count all disturbance pixel
    int n_prev = 0;
    int n_next = 0;
    
    patches.reset(); // clear list of patches
    
    for (int j = 0; j< dlength; ++j) {
      int dx = x + xdelta(j);
      int dy = y + ydelta(j);
      if ( dx>=0 && dx<ncol && dy>=0 && dy<nrow ) {
        // pixel values
        int pxval = r(dy,dx);
        if (pxval == value)
          ++n;
        if (pxval > 0)
          ++n_total;
        if (pxval == value - 1)
          ++n_prev;
        if (pxval == value + 1)
          ++n_next;
        // patches: flag each ID with a single bit
        int patch_val = patchdata(dy * ncol + dx);
        if (patch_val > 0 )
          patches.set(patch_val);
        
      } 
      
    }
    result(i,0) = n;
    result(i,1) = n_total;
    result(i,2) = n_prev;
    result(i,3) = n_next;
    result(i,4) = patches.count();
    
    //if (i % 1000 == 0) {
    //  Rcpp::checkUserInterrupt();
    //  Rcout << i << "/" << xidx.length() << " patches" << std::endl;
    //}
  }); 
  Rcout << "Finsihed! " << std::endl;
  
  return result;
  
}



// [[Rcpp::export]]
Rcpp::IntegerMatrix countPxAll(const Rcpp::IntegerMatrix &r, 
                               int lookup_value,
                               Rcpp::NumericVector xdelta,
                               Rcpp::NumericVector ydelta) {
  
  Rcpp::IntegerMatrix result( r.nrow(), r.ncol() );
  
  int ncol = r.ncol();
  int nrow = r.nrow();
  int dlength = xdelta.length();
  int n=0;
  
  for (int y = 0; y<nrow; ++y) {
    for (int x = 0; x<ncol; ++x) {
      if (r(y,x) == lookup_value) {
        ++n;
        for (int j = 0; j< dlength; ++j) {
          int dx = x + xdelta(j);
          int dy = y + ydelta(j);
          if ( dx>=0 && dx<ncol && dy>=0 && dy<nrow && r(dy, dx) == lookup_value ) 
            result(dy, dx)++;
        }
      }
    }
    if (y % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      Rcout << y << "/" << nrow << " lines, " << n << " px processed (total)" << std::endl;
    }
  }
  Rcout << "Finsihed! " << n << " px processed (total)" << std::endl;
  
  return result;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


/* Use:
 
 yr <- 1990
 
 out2 <- out %>% filter(year==yr)
 
 restab <- countPxCentroid3( dist_map_matrix, 
 out2$col, out2$row,
 out2$year,
 dtab$ix, dtab$iy,
 patches[[as.character(yr)]]@data@values)
 
 restab <- cbind( out2, restab)
 

 */