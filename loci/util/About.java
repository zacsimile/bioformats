//
// About.java
//

/*
LOCI common classes for use with VisBio, Bio-Formats, 4D Data Browser, etc.
Copyright (C) 2006 Curtis Rueden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Library General Public License for more details.

You should have received a copy of the GNU Library General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

package loci.util;

import javax.swing.JOptionPane;

/** Displays a small information dialog about this package. */
public abstract class About {

  public static void main(String[] args) {
    JOptionPane.showMessageDialog(null,
      "LOCI Common Classes\n" +
      "Built @date@\n\n" +
      "The LOCI common classes are LOCI software written by\n" +
      "Curtis Rueden, for use with other LOCI software packages.",
      "LOCI Common Classes", JOptionPane.INFORMATION_MESSAGE);
    System.exit(0);
  }

}
