import { defineConfig } from 'vite'

export default defineConfig({
  base: './',
  server: {
    port: 3000,
    open: true
  },
  build: {
    target: 'es2015',
    outDir: 'dist'
  }
})