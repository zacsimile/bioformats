package loci.formats.in;

import java.io.IOException;

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

    public DCIMGReader() {
        super("Hamamatsu DCIMG", "dcimg");
        suffixSufficient = false;
        domains = new String[] {FormatTools.UNKNOWN_DOMAIN};
    }

    @Override
    public boolean isThisType(RandomAccessInputStream stream) throws IOException {
        String desc = stream.readString(5);
        if (!desc.equals("DCIMG")) return false;
        return true;
    }
    
    @Override
    public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h) 
        throws FormatException, IOException
    {
        FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);

        // DCIMG is stored column major
        in.seek(sessionOffset + dataOffset + byteFactor*y*getSizeX());
        for (int row=h-1; row>=0; row--) {
            in.skipBytes(byteFactor*x);
            in.read(buf, byteFactor*row*w, byteFactor*w);
            in.skipBytes(byteFactor*(getSizeX() - w - x));
        }

        return buf;
    }

    @Override
    protected void initFile(String id) throws FormatException, IOException {
        super.initFile(id);
        in = new RandomAccessInputStream(id);

        in.order(IS_LITTLE);  // big endian 

        // confirm this is DCIMG
        String desc = in.readString(5);
        if (!desc.equals("DCIMG")) throw new FormatException("Not a valid DCIMG file.");

        in.skipBytes(3);

        long version = in.readUnsignedInt();  // DCIMG version number
        if ((!(version == DCAM_VERSION_0)) && (!(version == DCAM_VERSION_1))) {
            throw new FormatException(String.format("Unknown DCIMG version number %d.", version));
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
        m.sizeZ = 1;
        m.imageCount = getSizeZ() * getSizeC();

        if (version == DCAM_VERSION_0) {
            parseDCAMVersion0Header();
        } else if (version == DCAM_VERSION_1) {
            parseDCAMVersion1Header();
        }

        if (pixelType == DCIMG_PIXELTYPE_MONO8) {
            m.pixelType = FormatTools.UINT8;
            byteFactor = 1;
        } else if (pixelType == DCIMG_PIXELTYPE_MONO16) {
            m.pixelType = FormatTools.UINT16;
            byteFactor = 2;
        }

        addGlobalMeta("Version", version);

        // The metadata store we're working with.
        MetadataStore store = makeFilterMetadata();
        MetadataTools.populatePixels(store, this);

    }

    private void parseDCAMVersion0Header() throws IOException {
        
    }

    private void parseDCAMVersion1Header() throws IOException {
        CoreMetadata m = core.get(0);

        in.seek(sessionOffset);
        long sessionLength = in.readUnsignedInt();
        in.skipBytes(20);
        long pseudoOffset = in.readUnsignedInt(); 
        in.skipBytes(32);  // unknown numbers 1, 144, 65537
        m.sizeT = (int) in.readUnsignedInt();
        pixelType = in.readUnsignedInt();
        long mystery1 = in.readUnsignedInt();
        // TODO: Is this casting always legal? I think not...
        m.sizeX = (int) in.readUnsignedInt();  // num columns (this is a column-major format)
        m.sizeY = (int) in.readUnsignedInt();  // num rows
        long bytesPerRow = in.readUnsignedInt();
        long bytesPerImage = in.readUnsignedInt();
        in.skipBytes(8);
        dataOffset = in.readUnsignedInt();
        in.skipBytes(16);
        long bytesPerFrame = in.readUnsignedInt();

    }
}
