package loci.formats.in;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.meta.MetadataStore;

/**
 * DCIMGReader reads Hamamatsu DCIMG files.
 * 
 * Follows spec in https://github.com/python-microscopy/python-microscopy/blob/master/PYME/IO/dcimg.py.
 */
public class DCIMGReader extends FormatReader {

  // -- Constants --

  private static final boolean IS_LITTLE = true;

  private static final long DCAM_VERSION_0 = 0x7;
  private static final long DCAM_VERSION_1 = 0x1000000;

  private static final long DCIMG_PIXELTYPE_NONE = 0x00000000; // defined, but I don't know what to do with this
  private static final long DCIMG_PIXELTYPE_MONO8 = 0x00000001;
  private static final long DCIMG_PIXELTYPE_MONO16 = 0x00000002;

  // -- Fields --

  private long sessionOffset;
  private long dataOffset;
  private long pixelType;
  private int byteFactor;

  private List<String> companionFiles = new ArrayList<String>();
  private String[] uniqueFiles;

  public DCIMGReader() 
  {
    super("Hamamatsu DCIMG", "dcimg");
    suffixSufficient = false;
    domains = new String[] {FormatTools.UNKNOWN_DOMAIN};
  }

  @Override
  public boolean isSingleFile(String id) 
    throws FormatException, IOException 
  {
    return false;
  }

  @Override
  public int fileGroupOption(String id) 
    throws FormatException, IOException 
  {
    return FormatTools.CAN_GROUP;
  }

  @Override
  public String[] getSeriesUsedFiles(boolean noPixels) 
  {
    FormatTools.assertId(currentId, true, 1);

    if (!isGroupFiles() || noPixels) return null;

    return uniqueFiles;

  }

  @Override
  public boolean isThisType(RandomAccessInputStream stream) 
    throws IOException 
  {
    String desc = stream.readString(5);
    if (!desc.equals("DCIMG")) return false;
    return true;
  }
  
  @Override
  public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h) 
    throws FormatException, IOException
  {
    FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);

    RandomAccessInputStream stream = new RandomAccessInputStream(uniqueFiles[no]);

    // DCIMG is stored column major
    stream.seek(sessionOffset + dataOffset + byteFactor*y*getSizeX());
    for (int row=h-1; row>=0; row--) {
      stream.skipBytes(byteFactor*x);
      stream.read(buf, byteFactor*row*w, byteFactor*w);
      stream.skipBytes(byteFactor*(getSizeX() - w - x));
    }

    stream.close();
    return buf;
  }

  @Override
  protected void initFile(String id)
    throws FormatException, IOException
  {
    super.initFile(id);
    in = new RandomAccessInputStream(id);

    in.order(IS_LITTLE);  // little endian 

    // confirm this is DCIMG
    if (!isThisType(in)) throw new FormatException("Not a valid DCIMG file.");

    in.seek(8);

    long version = in.readUnsignedInt();  // DCIMG version number
    if ((!(version == DCAM_VERSION_0)) && (!(version >= DCAM_VERSION_1))) {
      throw new FormatException(String.format("Unknown DCIMG version number %d.", version));
    }

    if (version > DCAM_VERSION_1) {
      LOGGER.warn(String.format("Your file is DCAM version %d, but only %d is guaranteed to work.", version, DCAM_VERSION_1));
    }

    in.skipBytes(20);

    long numSessions = in.readUnsignedInt();
    long numFrames = in.readUnsignedInt();
    sessionOffset = in.readUnsignedInt();
    in.skipBytes(4);
    long fileSize = in.readUnsignedInt();
    in.skipBytes(12);
    long fileSize2 = in.readUnsignedInt();
    if (fileSize != fileSize2) throw new FormatException("Improper header. File sizes do not match.");
    in.skipBytes(16);
    long mystery1 = in.readUnsignedInt();  // 1024 in all examples

    CoreMetadata m = core.get(0);

    m.dimensionOrder = "XYZCT";
    m.rgb = false;
    m.interleaved = false;
    m.littleEndian = IS_LITTLE;
    m.indexed = false;
    m.falseColor = false;
    m.metadataComplete = true;
    m.thumbnail = false;
    m.sizeC = 1;

    if (version == DCAM_VERSION_0) {
      parseDCAMVersion0Header(in);
    } else if (version == DCAM_VERSION_1) {
      parseDCAMVersion1Header(in);
    }

    if (pixelType == DCIMG_PIXELTYPE_MONO8) {
      m.pixelType = FormatTools.UINT8;
      byteFactor = 1;
    } else if (pixelType == DCIMG_PIXELTYPE_MONO16) {
      m.pixelType = FormatTools.UINT16;
      byteFactor = 2;
    }

    if (isGroupFiles()) {
      Location currentFile = new Location(id).getAbsoluteFile();
      Location directory = currentFile.getParentFile();
      scanDirectory(directory);
    } else {
      String file = new Location(id).getAbsolutePath();
      companionFiles.add(file);
    }

    uniqueFiles = companionFiles.toArray(new String[companionFiles.size()]);
    m.sizeZ = uniqueFiles.length;
    m.imageCount = getSizeZ() * getSizeC();

    LOGGER.debug("sizeX: {} sizeY: {} sizeZ: {} sizeC: {} sizeT: {}", m.sizeX, m.sizeY, m.sizeZ, m.sizeC, m.sizeT);

    addGlobalMeta("Version", version);

    // The metadata store we're working with.
    MetadataStore store = makeFilterMetadata();
    MetadataTools.populatePixels(store, this);

    in.close();

  }

  /* Logic copied from DicomReader */
  private void scanDirectory(Location dir)
    throws FormatException, IOException 
  {
    String[] files = dir.list(true);
    if (files == null) return;
    Arrays.sort(files);
    for (String f : files) {
      String file = new Location(dir, f).getAbsolutePath();
      LOGGER.debug("Checking file {}", file);
      addFileToList(file);
    }
  }

  private void addFileToList(String file)
    throws FormatException, IOException
  {
    String ext = getExtension(file);
    if (!ext.equals("dcimg")) {
      LOGGER.debug("File {} with extension {} failed extension check", file, ext);
      return;
    }
    RandomAccessInputStream stream = new RandomAccessInputStream(file);
    if (!isThisType(stream)) {
      LOGGER.debug("File {} is not DCIMG", file);
      stream.close();
      return;
    }
    // stream.order(IS_LITTLE);
    // now check width/height
    companionFiles.add(file);
    stream.close();
  }

  private String getExtension(String file)
  {
    String ext = "";
    int i = file.lastIndexOf(".");
    if (i > 0) {
      ext = file.substring(i+1);
    }
    return ext;
  }

  private void parseDCAMVersion0Header(RandomAccessInputStream stream)
    throws IOException 
  {
    CoreMetadata m = core.get(0);

    stream.seek(sessionOffset);
    long sessionLength = stream.readUnsignedInt();
    stream.skipBytes(4);
    long pseudoOffset = stream.readUnsignedInt(); 
    stream.skipBytes(20);  
    m.sizeT = (int) stream.readUnsignedInt();
    pixelType = stream.readUnsignedInt();
    long mystery1 = stream.readUnsignedInt();
    // TODO: This will fail for large numbers
    m.sizeX = (int) stream.readUnsignedInt();  // num columns (this is a column-major format)
    long bytesPerRow = stream.readUnsignedInt();
    m.sizeY = (int) stream.readUnsignedInt();  // num rows
    long bytesPerImage = stream.readUnsignedInt();
    stream.skipBytes(8);
    dataOffset = stream.readUnsignedInt();
    long offsetToFooter = stream.readUnsignedInt();  /// TODO: Deal with footer
  }

  private void parseDCAMVersion1Header(RandomAccessInputStream stream)
    throws IOException 
  {
    CoreMetadata m = core.get(0);

    stream.seek(sessionOffset);
    long sessionLength = stream.readUnsignedInt();
    stream.skipBytes(20);
    long pseudoOffset = stream.readUnsignedInt(); 
    stream.skipBytes(32);  // unknown numbers 1, 144, 65537
    m.sizeT = (int) stream.readUnsignedInt();
    pixelType = stream.readUnsignedInt();
    long mystery1 = stream.readUnsignedInt();
    // TODO: This will fail for large numbers
    m.sizeX = (int) stream.readUnsignedInt();  // num columns (this is a column-major format)
    m.sizeY = (int) stream.readUnsignedInt();  // num rows
    long bytesPerRow = stream.readUnsignedInt();
    long bytesPerImage = stream.readUnsignedInt();
    stream.skipBytes(8);
    dataOffset = stream.readUnsignedInt();
    stream.skipBytes(16);
    long bytesPerFrame = stream.readUnsignedInt();
  }
}
