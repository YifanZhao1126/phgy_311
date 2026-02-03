import { defineConfig } from 'vite'

export default defineConfig({
  base: './',
  server: {
    port: 3000,
    open: true
  },
  build: {
    target: 'es2015',
    outDir: 'docs',
    rollupOptions: {
      output: {
        assetFileNames: 'assets/[name]-[hash][extname]',
        chunkFileNames: 'assets/[name]-[hash].js',
        entryFileNames: 'assets/[name]-[hash].js'
      }
    }
  }
})