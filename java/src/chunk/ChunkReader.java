package chunk;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

/**
 * A simple chunk reader class, similar to Python's Chunk library. [That
 * library's skip()-to-end method should be private.]
 *<p/>
 * CAUTION: If a read method raises an IOException, the offset into the
 * underlying input stream will be unknown!
 */
public class ChunkReader {
    static final Charset ASCII = Charset.forName("US-ASCII");

    public final String chunkType;
    public final int chunkSize;

    private DataInputStream input;
    private boolean align;
    private int offset;
    private boolean closed;

    /**
     * Read all the chunks in inputStream into a List<ChunkWriter>. (This is not
     * to imply that the caller will write them. ChunkWriter is just a POJO to
     * hold a chunk's contents.)
     *<p/>
     * TODO(jerry): An Iterator would be smarter than a List.
     */
    public static List<ChunkWriter> readAll(
            InputStream inputStream, boolean align) {
        List<ChunkWriter> result = new ArrayList<>();

        while (true) {
            try {
                ChunkReader reader = new ChunkReader(inputStream, align);
                ChunkWriter chunk =
                        new ChunkWriter(reader.getName(), reader.read());
                result.add(chunk);
            } catch (EOFException e) {
                break;
            } catch (IOException e) {
                e.printStackTrace();
                break;
            }
        }

        return result;
    }

    /**
     * Open a ChunkReader to read the next chunk from the given input stream.
     * This supports network byte order (big-endian) chunkSize fields, like
     * EA-IFF 85. An IFF reader would also need a method to open a nested chunk
     * reader using e.g. a bounded FilterInputStream onto inputStream.
     *
     * @param inputStream the stream to read from
     * @param align whether to read past a chunk alignment pad byte. Pass true
     *              for IFF; false for simpler files.
     *
     * @throws EOFException when there isn't another chunk to read
     * @throws IOException on other I/O errors
     */
    public ChunkReader(InputStream inputStream, boolean align)
            throws IOException {
        byte[] chunkTypeBytes = new byte[4];

        this.input = new DataInputStream(inputStream);
        this.align = align;
        this.input.readFully(chunkTypeBytes);
        this.chunkType = new String(chunkTypeBytes, ASCII);
        this.chunkSize = this.input.readInt();
        this.offset = 0;
        this.closed = false;
    }

    public String getName() {
        return chunkType;
    }

    public int getSize() {
        return chunkSize;
    }

    /**
     * Close the ChunkReader, seeking the input stream to the next chunk (past
     * the rest of this chunk body and any alignment pad byte) so the caller can
     * open the next ChunkReader. Idempotent.
     */
    public void close() throws IOException {
        if (!closed) {
            int pad = align && (chunkSize & 1) == 1 ? 1 : 0;

            input.skipBytes(chunkSize - offset + pad);
            offset = chunkSize;
            closed = true;
            input = null;
        }
    }

    /**
     * Seek in the chunk body. See {@link #skipBytes(int)}
     *
     * TODO(jerry): Can input.skipBytes() seek backwards?
     *
     * @param pos number of bytes relative to the...
     * @param whence 0 = chunk start, 1 = current offset, 2 = chunk end.
     * @return the number of bytes skipped over.
     */
    public int seek(int pos, int whence) throws IOException {
        int n;

        if (closed) {
            throw new IOException("Can't seek in a closed chunk");
        }

        if (whence == 0) {
            n = pos - offset;
        } else if (whence == 2) {
            n = chunkSize + pos - offset;
        } else {
            n = pos;
        }

        if (offset + n < 0 || offset + n > chunkSize) {
            throw new IOException("Can't seek out of chunk range");
        }

        int actual = input.skipBytes(n);
        offset += actual;
        return actual;
    }

    /** Return the current read-offset into the chunk body. */
    public int tell() {
        return offset;
    }

    /** Read the remainder of the chunk body. */
    public byte[] read() throws IOException {
        return read(chunkSize - offset);
    }

    /**
     * Read the given number of bytes from the chunk body, capped to the number
     * of bytes remaining.
     */
    public byte[] read(int length) throws IOException {
        if (length > chunkSize - offset) {
            length = chunkSize - offset;
        }

        byte[] buffer = new byte[length];

        readFully(buffer);
        return buffer;
    }

    // --- Part of the DataInput interface ---

    /** Read b.length bytes into b[0:]. */
    public void readFully(byte[] b) throws IOException {
        readFully(b, 0, b.length);
    }

    /** Read len bytes into b[off:off+len]. */
    public void readFully(byte[] b, int off, int len) throws IOException {
        if (closed) {
            throw new IOException("Can't read a closed chunk");
        }

        if (len < 0 || offset + len > chunkSize) {
            throw new EOFException("read length is out of bounds");
        }

        input.readFully(b, off, len);
        offset += len;
    }

    /** Skip over some bytes in the chunk body. */
    public int skipBytes(int n) throws IOException {
        return seek(n, 1);
    }
}
