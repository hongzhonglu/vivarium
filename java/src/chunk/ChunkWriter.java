package chunk;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/** A simple chunk writer. */
public class ChunkWriter {
    public final String chunkType;
    public final byte[] body;

    public ChunkWriter(String chunkType, byte[] body) {
        this.chunkType = chunkType;
        this.body = body;
    }

    public int size() { return body.length; }

    /** Write the chunk to the output stream. */
    public void write(OutputStream outputStream, boolean align) throws IOException {
        DataOutputStream output = new DataOutputStream(outputStream);
        byte[] chunkTypeBytes = chunkType.getBytes(ChunkReader.ASCII);

        output.write(chunkTypeBytes);
        output.writeInt(body.length);
        output.write(body);
        if (align && (body.length & 1) == 1) {
            output.writeByte(0);
        }
        output.flush();
    }
}
